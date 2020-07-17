import copy
import subprocess
import statistics
import operator
import os
import logging

from typing import Optional, Dict, List

from .sglib import (
    StarAllele,
    SNPAllele,
    VCFFile,
    read_gene_table,
    read_snp_table,
    build_snpdb,
    read_star_table,
    build_stardb,
    read_phenotype_table,
    parse_vcf_fields,
    vcf2biosamples,
    sort_star_names,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def cyp2a6(sample, stardb):
    call_sv1(sample, "gc_e1e4", "*34", stardb)
    call_sv1(sample, "gc_e1e2", "*12", stardb)
    call_tandem(sample, "dup2", "*2", "*S1", stardb)
    call_tandem(sample, "dup2", "*1", "*S1", stardb)
    call_tandem(sample, "gc_e9", "*1", "*S2", stardb)
    call_tandem(sample, "dup7", "*1", "*S3", stardb)
    call_tandem(sample, "dup7b", "*1", "*S6", stardb)
    cyp2a6_svcomb(sample, stardb)

def cyp2b6(sample, stardb):
    call_sv1(sample, "gc_i4e9", "*29", stardb)

def cyp2d6(sample, stardb):
    call_tandem(sample, 'del1', '*S2', '*1', stardb, ordered = True)
    call_tandem(sample, 'gc_i1e9_5', '*68x5', '*4', stardb)
    call_tandem(sample, 'gc_i1e9', '*68', '*4', stardb)
    call_tandem(sample, 'gc_i1e9', '*S1', '*1', stardb)
    call_tandem(sample, 'gc_e9', '*4N', '*4', stardb)
    call_tandem(sample, 'gc_e9_3', '*36x3', '*10', stardb)
    call_tandem(sample, 'gc_e9_7', '*36x7', '*10', stardb)
    call_tandem(sample, 'gc_e9', '*36', '*10', stardb)
    call_tandem(sample, 'gc_e9', '*83', '*2', stardb)
    call_tandem(sample, 'gc_7to6_i4', '*13A', '*2', stardb)
    call_sv1(sample, 'gc_e1e7', '*13C', stardb)
    call_sv1(sample, 'gc_7to6_i1', '*13B', stardb)
    cyp2d6_svcomb(sample, stardb)

def cyp2e1(sample, stardb):
    call_sv1(sample, "dup_e7e9", "*S1", stardb)

def gstm1(sample, stardb):
    gstm1_svcomb(sample, stardb)

def gstt1(sample, stardb):
    gstt1_svcomb(sample, stardb)

def slc22a2(sample, stardb):
    call_sv1(sample, "del_i9", "*S1", stardb)
    call_sv1(sample, "del_e11", "*S2", stardb)

def slco1b1(sample, stardb):
    call_sv1(sample, "dup1", "*S3", stardb)

def ugt1a4(sample, stardb):
    call_sv1(sample, 'del_i1', '*S1', stardb)
    call_sv1(sample, 'del2', '*S2', stardb)

def ugt2b15(sample, stardb):
    call_sv1(sample, "del_i3e6", "*S1", stardb)

def ugt2b17(sample, stardb):
    ugt2b17_svcomb(sample, stardb)

def new_tandem(sv, star_names, stardb):
    stars = [stardb[x] for x in star_names]
    core = list(set([x for y in [z.core for z in stars] for x in y]))
    tandem = StarAllele()
    tandem.name, tandem.score, tandem.core, tandem.sv = '+'.join(star_names), sum([x.score for x in stars]), copy.deepcopy(core), sv
    return tandem
    
def new_dup(sX, cnv):
    times = int(cnv.replace('cnv', ''))
    dup = StarAllele()
    dup.name, dup.score, dup.core, dup.sv = sX.name + 'x' + str(times), sX.score * times, copy.deepcopy(sX.core), cnv
    return dup

def remove_select(hap, stars):
    '''
    Args:
        hap (list of Stars)
        stars (list of Stars)
    '''
    for i in reversed(range(len(hap))):
        if hap[i].name in [x.name for x in stars]:
            del hap[i]

def remove_sv(hap, l = []):
    for i in reversed(range(len(hap))):
        if hap[i].sv:
            if hap[i].name in l:
                continue
            else:
                del hap[i]
            
def which_has(sample, stars):
    '''
    Args:
        sample (Sample)
        stars (list of str)
    Returns:
        i (int)
    '''
    h1 = set(stars).issubset([x.name for x in sample.hap[0].cand])
    h2 = set(stars).issubset([x.name for x in sample.hap[1].cand])
    if h1 and h2:
        i = 3
    elif h1 and not h2:
        i = 1
    elif not h1 and h2:
        i = 2
    else:
        i = 0
    return i

def call_sv1(sample, sv, x, stardb):
    '''
    This function calls a sample's final genotype if the sample has only one SV.
    If a SV-carrying allele is in LD with other alleles, the funtion takes those other alleles as input. 
    x = the name of the star allele with the SV
    '''
    if sample.gt or sample.sv != ["no_sv", sv]:
        return
    if stardb[x].core:
        h1 = set(stardb[x].core).issubset(sample.hap[0].obs)
        h2 = set(stardb[x].core).issubset(sample.hap[1].obs)
    else:
        h1 = True
        h2 = True
    if not h1 and not h2:
        return
    elif h1 and not h2:
        i, j = 1, 0
    elif not h1 and h2:
        i, j = 0, 1
    else:        
        l1 = copy.deepcopy(sample.hap[0].cand)
        l2 = copy.deepcopy(sample.hap[1].cand)
        remove_sv(l1)
        remove_sv(l2)
        if l1[0].name == l2[0].name:
            i, j = 0, 1
        else:
            l3 = sorted([l1[0], l2[0]], key = lambda x: x.rank)
            if l3[0].name == l1[0].name:
                i, j = 0, 1
            else:
                i, j = 1, 0
    remove_sv(sample.hap[i].cand)
    remove_sv(sample.hap[j].cand, [x])
    sample.gt = True

def call_tandem(sample, sv, x, y, stardb, ordered = False):
    """
    Calls a tandem duplication allele containing two gene copies (e.g., CYP2D6*36+*10)
    x = the name of the 1st star allele in the tandem (e.g., '*36')
    y = the name of the 2nd star allele in the tandem (e.g., '*10')
    """
    if sample.gt or sample.sv != ["no_sv", sv]:
        return        
    h1 = set([x, y]).issubset([_.name for _ in sample.hap[0].cand])
    h2 = set([x, y]).issubset([_.name for _ in sample.hap[1].cand])
    if not h1 and not h2:
        return
    elif h1 and not h2:
        i, j = 0, 1
    elif not h1 and h2:
        i, j = 1, 0
    else:
        l1 = copy.deepcopy(sample.hap[0].cand)
        l2 = copy.deepcopy(sample.hap[1].cand)
        remove_sv(l1)
        remove_sv(l2)
        if l1[0] == l2[0]:
            i, j = 1, 0
        else:
            l3 = sorted([l1[0], l2[0]], key = lambda x: x.rank)
            if l3[0] == l1[0]:
                i, j = 1, 0    
            else:
                i, j = 0, 1
    sX = stardb[x]
    sY = stardb[y]
    
    # find SNPs shared by both star alleles
    overlap = []
    for snp in sX.core:
        if snp in sY.core:
            overlap.append(snp)
    
    # return if allele fraction in any of the shared SNPs is less than 0.4
    for snp in overlap:
        if [x for x in sample.hap[i].obs if x == snp][0].af < 0.4:
            return    
    
    tandem = new_tandem(sv, [sX.name, sY.name], stardb)
    sample.hap[i].cand.insert(0, tandem)
    remove_select(sample.hap[i].cand, [sX, sY])
    remove_sv(sample.hap[i].cand, [tandem.name])
    remove_sv(sample.hap[j].cand)
    if ordered:
        sample.hap[i].cand.sort(key = lambda x: x.rank)
    sample.gt = True

def call_cnv3(sample):
    if sample.gt or sample.sv != ["no_sv", 'cnv2']:
        return
    sX = sample.hap[0].cand[0]
    sY = sample.hap[1].cand[0]

    # Simplest case (e.g., *2/*2x2)
    if sX.name == sY.name:
        sample.hap[0].add_dup(2)
        sample.gt = True
        return

    sX_gene = sample.hap[0].af_mean_gene
    sY_gene = sample.hap[1].af_mean_gene
    if sX_gene == -1:
        sX_gene = 1 - sY_gene
    if sY_gene == -1:
        sY_gene = 1 - sX_gene
    diff_gene = sY_gene - sX_gene

    sX_main = sample.hap[0].af_mean_main
    sY_main = sample.hap[1].af_mean_main
    if sX_main == -1:
         sX_main = 1 - sY_main
    if sY_main == -1:
        sY_main = 1 - sX_main
    diff_main = sY_main - sX_main    
    
    fit_maf1, _ = fit_data(3, sample.hap[0].af_mean_gene)
    fit_maf2, _ = fit_data(3, sample.hap[1].af_mean_gene)
    
    means = [round(fit_maf1, 2), round(fit_maf2, 2)]

    f = lambda a, b: ((a == b) & (a == 0)) | (a * b > 0)
        
    if f(diff_gene, diff_main):
        if means == [0.33, 0.67]:
            sample.hap[1].add_dup(2)
            sample.gt = True
        elif means == [0.67, 0.33]:
            sample.hap[0].add_dup(2)
            sample.gt = True
    else:
        
        if abs(diff_main) > abs(diff_gene):
            
            if sY_main > sX_main:
                sample.hap[1].add_dup(2)
                sample.gt = True
                
            else:
                sample.hap[0].add_dup(2)
                sample.gt = True
        
        else:
            if means == [0.33, 0.67]:
                sample.hap[1].add_dup(2)
                sample.gt = True
            elif means == [0.67, 0.33]:
                sample.hap[0].add_dup(2)
                sample.gt = True

def fit_data(total_cn, af_mean_gene):
    """Return the fit MAF and CN."""
    maf_choices = []
    for i in range(1, total_cn):
        maf_choices.append(i / total_cn)
    fit_maf = min(maf_choices, key = lambda x: abs(x - af_mean_gene))
    fit_cn = maf_choices.index(fit_maf) + 1
    return fit_maf, fit_cn

def call_cnv_plus(sample):
    '''This function calls a final genotype with CN > 3 gene copies.'''
    if sample.gt:
        return
    if 'cnv' not in sample.sv[1]:
        return
    if sample.sv[0] == 'no_sv' and (sample.sv[1] == 'cnv0' or sample.sv[1] == 'cnv2'):
        return
    if sample.sv[0] == 'no_sv' and 'cnv' in sample.sv[1]:
        total_cn = int(sample.sv[1].replace('cnv', '')) + 1
    elif 'cnv' in sample.sv[0] and 'cnv' in sample.sv[1]:
        total_cn = int(sample.sv[0].replace('cnv', '')) + int(sample.sv[1].replace('cnv', ''))
        if total_cn < 4:
            return

    # allele fraction profile is not informative -- i.e. it's empty
    if sample.hap[0].af_mean_gene == -1:
        sample.hap[0].add_dup(total_cn - 1)
        sample.gt = True
        return

    fit_maf, fit_cn = fit_data(total_cn, sample.hap[0].af_mean_gene)
    sample.hap[0].add_dup(fit_cn)
    sample.hap[1].add_dup(total_cn - fit_cn)
    sample.gt = True

def    cyp2a6_svcomb(sample, stardb):
    if sample.gt:
        return
    gt = []
    for sv in sample.sv:
        if sv == 'cnv0':
            gt.append(stardb['*4'])
        elif sv == 'cnv2':
            gt.append(new_dup(stardb['*1'], sv))
        elif sv == 'gc_e1e2':
            gt.append(stardb['*12'])
        elif sv == 'gc_e1e4':
            gt.append(stardb['*34'])
        elif sv == 'dup2':
            gt.append(new_tandem(sv, ['*1', '*S1'], stardb))
        elif sv == 'gc_e9':
            gt.append(new_tandem(sv, ['*1', '*S2'], stardb))
        elif sv == 'dup7':
            gt.append(new_tandem(sv, ['*1', '*S3'], stardb))
        elif sv == 'dup7x2':
            gt.append(new_tandem(sv, ['*1', '*S3', '*S3'], stardb))
        elif sv == 'dup7b':
            gt.append(new_tandem(sv, ['*1', '*S6'], stardb))
    if len(gt) == 2:
        sample.hap[0].cand = [gt[0]]
        sample.hap[1].cand = [gt[1]]
        sample.gt = True    

def    cyp2d6_svcomb(sample, stardb):
    if sample.gt:
        return
    gt = []
    for sv in sample.sv:
        if sv == 'cnv0':
            gt.append(stardb['*5'])
        elif sv == 'gc_i1e9' and which_has(sample, ['*68', '*4']):
            gt.append(new_tandem(sv, ['*68', '*4'], stardb))
        elif sv == 'gc_i1e9' and which_has(sample, ['*S1', '*1']):
            gt.append(new_tandem(sv, ['*S1', '*1'], stardb))
        elif sv == 'gc_e9' and which_has(sample, ['*4N', '*4']):
            gt.append(new_tandem(sv, ['*4N', '*4'], stardb))
        elif sv == 'gc_e9' and which_has(sample, ['*36', '*10']):
            gt.append(new_tandem(sv, ['*36', '*10'], stardb))
        elif sv == 'gc_e9' and which_has(sample, ['*83', '*2']):
            gt.append(new_tandem(sv, ['*83', '*2'], stardb))
        elif sv == 'gc_7to6_i4' and which_has(sample, ['13A', '*2']): 
            gt.append(new_tandem(sv, ['13A', '*2'], stardb))
        elif sv == 'gc_7to6_i1':
            gt.append(stardb['*13B'])
        elif sv == 'gc_e1e7':
            gt.append(stardb['*13C'])
    cnv = None
    for sv in sample.sv:
        if 'cnv' in sv and sv != 'cnv0':
            cnv = sv
    if cnv:
        if '*68+*4' in [x.name for x in gt]:
            svcomb_tandem_cnv(gt, sample, ['*68', '*4'], cnv)
        elif '*S1+*1' in [x.name for x in gt]:
            svcomb_tandem_cnv(gt, sample, ['*S1', '*1'], cnv)
        elif '*4N+*4' in [x.name for x in gt]:
            svcomb_tandem_cnv(gt, sample, ['*4N', '*4'], cnv)
        elif '*36+*10' in [x.name for x in gt]:
            svcomb_tandem_cnv(gt, sample, ['*36', '*10'], cnv)
        elif '*83+*2' in [x.name for x in gt]:
            svcomb_tandem_cnv(gt, sample, ['*83', '*2'], cnv)
        elif '*13A+*2' in [x.name for x in gt]:
            svcomb_tandem_cnv(gt, sample, ['*13A', '*2'], cnv)
        elif '*13B' in [x.name for x in gt]:
            svcomb_sv1_cnv(gt, sample, '*13B', cnv)

    if len(gt) == 2:
        sample.hap[0].cand = [gt[0]]
        sample.hap[1].cand = [gt[1]]
        sample.gt = True

def svcomb_sv1_cnv(gt, sample, sX_name, cnv):
    i = which_has(sample, [sX_name])
    if not i:
        return
    if i != 3:
        j = {0: 1, 1: 0}[i - 1]
    elif i == 3:
        j = 0
    l = copy.deepcopy(sample.hap[j].cand)
    remove_sv(l)
    sY = new_dup(l[0], cnv)
    gt.insert(j, sY)

def svcomb_tandem_cnv(gt, sample, tandem, cnv):
    i = which_has(sample, [tandem[0], tandem[1]])
    if not i:
        return
    if i != 3:
        j = {0: 1, 1: 0}[i - 1]
        l = copy.deepcopy(sample.hap[j].cand)
        remove_sv(l)
        sX = new_dup(l[0], cnv)
        gt.insert(j, sX)
    elif i == 3:
        for x in sample.hap[0].cand:
            if x.name == tandem[1]:
                gt.append(new_dup(x, cnv))
                break

def    gstt1_svcomb(sample, stardb):
    if sample.gt:
        return
    gt = []
    for sv in sample.sv:
        if sv == 'cnv0':
            gt.append(stardb["*2"])
    if len(gt) == 2:
        sample.hap[0].cand = [gt[0]]
        sample.hap[1].cand = [gt[1]]
        sample.gt = True

def    gstm1_svcomb(sample, stardb):
    if sample.gt:
        return
    gt = []
    for sv in sample.sv:
        if sv == 'cnv0':
            gt.append(stardb["*2"])
    if len(gt) == 2:
        sample.hap[0].cand = [gt[0]]
        sample.hap[1].cand = [gt[1]]
        sample.gt = True    

def    ugt2b17_svcomb(sample, stardb):
    if sample.gt:
        return
    gt = []
    for sv in sample.sv:
        if sv == 'cnv0':
            gt.append(stardb["*2"])
    if len(gt) == 2:
        sample.hap[0].cand = [gt[0]]
        sample.hap[1].cand = [gt[1]]
        sample.gt = True

def remove_extra_s1(biosamples):
    def f(l):
        if len(l) == 1:
            return
        for i in reversed(range(len(l))):
            if l[i].name == '*1':
                del l[i]

    for biosample in biosamples:
        f(biosample.hap[0].cand)
        f(biosample.hap[1].cand)
        f(biosample.dip_cand)

def write_result_file(biosamples, out):
    list2str = lambda x: '.' if not x else ','.join([str(x) for x in x])

    header = [
        'name', 'status', 'hap1_main', 'hap2_main', 'hap1_cand', 'hap2_cand',
        'hap1_score', 'hap2_score', 'dip_score', 'phenotype', 'dip_sv',
        'hap1_sv', 'hap2_sv', 'ssr', 'dip_cand', 'hap1_main_core',
        'hap2_main_core', 'hap1_main_tag', 'hap2_main_tag',
        'hap1_af_mean_gene', 'hap2_af_mean_gene', 'hap1_af_mean_main',
        'hap2_af_mean_main'
    ]

    with open(f"{out}/genotype.txt", "w") as f:
        f.write('\t'.join(header) + '\n')
        for biosample in biosamples:
            fields = ['.' for x in header]
            status = 'g' if biosample.gt else 'ng'
            fields[header.index('name')] = biosample.name
            fields[header.index('status')] = status
            if status == 'g': fields[header.index('hap1_main')] = biosample.hap[0].cand[0].name
            if status == 'g': fields[header.index('hap2_main')] = biosample.hap[1].cand[0].name
            fields[header.index('hap1_cand')] = ','.join([x.name for x in biosample.hap[0].cand])
            fields[header.index('hap2_cand')] = ','.join([x.name for x in biosample.hap[1].cand])
            if status == 'g': fields[header.index('hap1_score')] = biosample.hap[0].cand[0].score
            if status == 'g': fields[header.index('hap2_score')] = biosample.hap[1].cand[0].score
            if status == 'g': fields[header.index('dip_score')] = biosample.hap[0].cand[0].score + biosample.hap[1].cand[0].score
            if status == 'g': fields[header.index('phenotype')] = biosample.pt
            fields[header.index('dip_sv')] = ','.join(biosample.sv)
            if status == 'g': fields[header.index('hap1_sv')] = biosample.hap[0].sv
            if status == 'g': fields[header.index('hap2_sv')] = biosample.hap[1].sv
            fields[header.index('ssr')] = biosample.ssr
            fields[header.index('dip_cand')] = ','.join([x.name for x in biosample.dip_cand])
            fields[header.index('hap1_main_core')] = list2str([x.summary() for x in biosample.hap[0].obs if x in biosample.hap[0].cand[0].core])
            fields[header.index('hap2_main_core')] = list2str([x.summary() for x in biosample.hap[1].obs if x in biosample.hap[1].cand[0].core])
            fields[header.index('hap1_main_tag')] = list2str([x.summary() for x in biosample.hap[0].obs if x in biosample.hap[0].cand[0].tag])
            fields[header.index('hap2_main_tag')] = list2str([x.summary() for x in biosample.hap[1].obs if x in biosample.hap[1].cand[0].tag])
            fields[header.index('hap1_af_mean_gene')] = biosample.hap[0].af_mean_gene
            fields[header.index('hap2_af_mean_gene')] = biosample.hap[1].af_mean_gene
            fields[header.index('hap1_af_mean_main')] = biosample.hap[0].af_mean_main
            fields[header.index('hap2_af_mean_main')] = biosample.hap[1].af_mean_main
            f.write('\t'.join([str(x) for x in fields]) + '\n')

def predict_phenotypes(biosamples, tg, pd):
    operator_dict = {
        "<": operator.lt,
        "<=": operator.le,
        ">": operator.gt,
        ">=": operator.ge,
        "==": operator.eq,
    }

    phenotype_table = read_phenotype_table(f"{pd}/phenotype_table.txt")

    for biosample in biosamples:
        for pt in phenotype_table[tg]:
            rules = phenotype_table[tg][pt]["rules"].strip(",").split(",")
            
            found = True
            for rule in rules:
                op = rule.split(":")[0]

                a = biosample.hap[0].cand[0].score
                b = biosample.hap[1].cand[0].score
                rule_score = float(rule.split(":")[1])

                if not operator_dict[op](a + b, rule_score):
                    found = False
                    break

            if found:
                biosample.pt = pt
                break

def assess_vcf(input_vcf):
    n_tt = len(input_vcf.data)
    n_ad = 0 # with AD field
    n_fp = 0 # fully phased
    n_np = 0 # not phased
    n_pp = 0 # partially phased

    for fields in input_vcf.data:
        v = parse_vcf_fields(fields)

        # Raise an error if the GT field is not found in the FORMAT column.
        if 'GT' not in v['format']:
            raise ValueError('GT field not found: ' + v['pos'])

        # Count the number of markers with the AD field.
        if 'AD' in v['format']:
            n_ad += 1

        # This function returns GT field separator for each sample.
        def f(x):
            gt = x.split(':')[0]
            if '/' in gt:
                return '/'
            elif '|' in gt:
                return '|'
            elif v['chrom'] in ['X', 'Y']:
                return
            else:
                raise ValueError('GT field separator not found: ' + v['pos'])

        # Determine the status of GT field separator for the marker.
        seps = set([f(x) for x in fields[9:] if f(x)])
        if len(seps) == 1 and '|' in seps:
            n_fp += 1
        elif len(seps) == 1 and '/' in seps:
            n_np += 1
        elif len(seps) == 2 and '|' in seps and '/' in seps:
            n_pp += 1
        else:
            raise ValueError('GT field not recognized: %s' % v['pos'])

    # Determine whether the AD field will be used.
    if n_tt > 0 and n_ad / n_tt > 0.8:
        vcf_ad = True
    else:
        vcf_ad = False

    # Get the GT field separator.
    if n_tt == 0:
        vcf_sep = None
    elif n_fp == n_tt:
        vcf_sep = '|'
    elif n_np == n_tt:
        vcf_sep = '/'
    else:
        vcf_sep = 'b'

    return vcf_sep, vcf_ad

def process_vcf(input_vcf, vcf_sep, vcf_ad, tg, gb):
    processed_vcf = input_vcf.copy(["header"])

    processed_vcf.meta = [
        f'##target_gene={tg}\n',
        f'##genome_build={gb}\n',
        '##fileformat=VCFv4.2\n',
        '##INFO=<ID=PS,Number=1,Type=String,Description="Phasing Status '
            '(A, in preparation; B1, ready for phasing as is; '
            'B2, ready for phasing after conformation to reference VCF; '
            'C1, excluded from phasing because marker is absent in '
            'reference VCF; C2, excluded from phasing because marker '
            'has different REF allele; C3, excluded from phasing because '
            'marker has no overlapping ALT alleles; D1, statistically '
            'phased; D2, manually phased with certainty; D3, already phased; '
            'E, omitted during statistical phasing)">\n',
        '##INFO=<ID=FE,Number=A,Type=String,Description='
            f'"Functional Effect">\n',
        '##INFO=<ID=RV,Number=A,Type=String,Description='
            f'"Reverting Variation">\n',
        '##INFO=<ID=SO,Number=A,Type=String,Description='
            f'"=Sequence Ontology">\n',
        '##INFO=<ID=VI,Number=A,Type=String,Description='
            f'"Variant Impact">\n',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths '
            'for the ref and alt alleles in the order listed">\n',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        '##FORMAT=<ID=HE,Number=4,Type=Integer,Description="Matching scores '
            'computed by the phase-by-extension algorithm">\n'
    ]

    for fields in input_vcf.data:
        v = parse_vcf_fields(fields)

        # Filter problematic markers.
        if v['ref'] == 'I' or v['ref'] == 'D' or '.' in v['alt'] or 'D' in v['alt']:
            logger.warning('Marker filtered due to invalid allele: %s' % v['pos'])
            continue
        
        missingness = ['.' in x.split(':')[0] for x in fields[9:]].count(True) / len(fields[9:])
        if missingness > 0.5:
            logger.warning('Marker filtered due to high missingness: {} ({:.2f})'.format(v['pos'], missingness))
            continue

        # Update the QUAL and FILTER data.
        fields[5] = '.'
        fields[6] = '.'

        # Update the INFO data.
        if vcf_sep == '|':
            fields[7] = 'PS=D3'
        else:
            fields[7] = 'PS=A'

        # Update the FORMAT data.
        if vcf_ad:
            fields[8] = 'GT:AD'
        else:
            fields[8] = 'GT'

        # Get the AD field index.
        if 'AD' in v['format']:
            i_ad = v['format'].index('AD')
        else:
            i_ad = None

        # Get the DP field index.
        if 'DP' in v['format']:
            i_dp = v['format'].index('DP')
        else:
            i_dp = None

        # This function returns the GT and AD fields.
        def f(x):
            gt = x.split(':')[0]

            # Determine the GT field separator.
            if '/' in gt:
                sep = '/'
            elif '|' in gt:
                sep = '|'
            else:
                sep = None

            # Remove phasing if the input VCF file is partially phased.
            if vcf_sep == 'b' and sep == '|':
                if gt == '.|.':
                    gt = './.'
                else:
                    gt = '/'.join(sorted(gt.split('|'), key=lambda x: int(x)))

            # Conform the GT field for the sex chromosomes.
            if not sep and v['chrom'] in ['X', 'Y']:
                if vcf_sep == '|':
                    if gt == '.':
                        gt = '.|.'
                    else:
                        gt = f'0|{gt}'
                else:
                    if gt == '.':
                        gt = './.'
                    else:
                        gt = f'0/{gt}'

            # Add the AD field.
            if vcf_ad and (gt == './.' or not i_ad):
                ad = ','.join(['0'] * (len(v['alt']) + 1))
            elif vcf_ad and i_ad:
                ad = x.split(':')[i_ad]
                # Variable-length AD field from the Genalice program (e.g. GT:AD:DP -> 1/1:21:22).
                if len(ad.split(',')) == 1 and sep and gt.split(sep)[0] == gt.split(sep)[1] and gt.split(sep)[0] != '0' and len(v['alt']) == 1:
                    if i_dp:
                        dp = x.split(':')[i_dp]
                        ad = str(int(dp) - int(ad)) + ',' + ad
                    else:
                        ad = f'0,{ad}'
                # Invalid AD field from the bcftools program (e.g. GT:AD -> 0/2:2,16).
                if len(ad.split(',')) != len(v['alt']) + 1:
                    ad = ','.join(['0'] * (len(v['alt']) + 1))
            else:
                ad = None

            # Return the genotype fields.
            if ad:
                return f'{gt}:{ad}'
            else:
                return gt

        fields[9:] = [f(x) for x in fields[9:]]
        processed_vcf.data.append(fields)

    return processed_vcf

def adjust_vcf(processed_vcf, stardb):
    """
    Conform multiallelic loci containing more than one indels to the star 
    allele table. For example, the UGT1A1 gene has a multiallelic locus 
    which defines three star alleles:

        *28 (234668879:CAT>CATAT)
        *36 (234668879:CAT>C)
        *37 (234668879:CAT>CATATAT)

    This locus can be represented in a VCF file in many different ways:

    VCF file					POS			REF		ALT
    1st VCF with *28 only		234668879	C		CAT
    2nd VCF with *36 only		234668879	CAT		C
    3rd VCF with *37 only		234668879	C		CATAT
    4th VCF with *28, *36		234668879	CAT		CATAT,C
    5th VCF with *28, *37		234668879	C		CAT,CATAT
    6th VCF with *28, *36, *37	234668879	CAT		CATAT,C,CATATAT
    ---------------------------------------------------------------
    Ref VCF has *28 only		234668879	C		CAT

    All VCF records will be adjusted to 234668879:CAT>CATAT,C,CATATAT.
    """

    adjusted_vcf = processed_vcf.copy(["meta", "header", "data"])

    for i in range(len(adjusted_vcf.data)):
        fields = adjusted_vcf.data[i]
        v = parse_vcf_fields(fields)
        
        # Skip if the locus does not have any indels.
        if len(v['ref']) == 1 and all([len(x) == 1 for x in v['alt']]):
            continue

        # Skip if the locus does not define any star alleles.
        star_list = []
        for name, star in stardb.items():
            if v['pos'] in [x.pos for x in star.core]:
                star_list.append(star)
        if not star_list:
            continue

        # Skip if the locus is already formatted properly.
        bool_list = []
        for x in v['alt']:
            is_found = False
            for star in star_list:
                if f"{v['pos']}:{v['ref']}>{x}" in [
                        f'{x.pos}:{x.hg}>{x.var}' for x in star.core
                    ]:
                    is_found = True
                    break
            bool_list.append(is_found)
        if all(bool_list):
            continue

        # Find matching alleles for the REF and ALT fields.
        snp_list = []
        for star in star_list:
            snp_list += star.core
        snp_list = list(set(snp_list))
        new_ref = []
        new_alt = []
        for x in v['alt']:
            result = [y for y in snp_list 
                if y.pos == v['pos'] and 
                (len(y.hg) - len(v['ref']) == len(y.var) - len(x))]
            if result:
                new_ref.append(result[0].hg)
                new_alt.append(result[0].var)
        new_ref = list(set(new_ref))

        # Skip if there are more than one allele for the REF field.
        if len(new_ref) > 1:
            continue

        # Skip if the number of alleles does not match for the ALT field.
        if len(new_alt) != len(v['alt']):
            continue

        # Update the REF and ALT fields.
        fields[3] = new_ref[0]
        fields[4] = ','.join(new_alt)
        logger.warning(f'Multiallelic marker reformatted: %s' % v['pos'])

    return adjusted_vcf

def conform_vcf(adjusted_vcf, ref_vcf, vcf_ad):    
    conformed_vcf = adjusted_vcf.copy(["meta", "header"])

    for fields1 in adjusted_vcf.data:
        v1 = parse_vcf_fields(fields1)
        is_found = False

        for i in range(len(ref_vcf.data)):
            fields2 = ref_vcf.data[i]
            v2 = parse_vcf_fields(fields2)

            # Skip if CHROM or POS does not match.
            if v1['chrom'] != v2['chrom'] or v1['pos'] != v2['pos']:
                continue

            # If the only difference is REF, check the next line just in case.
            if v1['ref'] != v2['ref']:
                fields3 = ref_vcf.data[i + 1]
                v3 = parse_vcf_fields(fields3)
                if v1['pos'] == v3['pos'] and v1['ref'] == v3['ref']:
                    continue
                else:
                    fields1[7] = 'PS=C2'
                    is_found = True
                    break

            # Skip because there are no matching ALT alleles.
            if not set(v1['alt']) & set(v2['alt']):
                is_found = True
                fields1[7] = 'PS=C3'
                break

            # Skip because there is no need for conformation.
            if v1['alt'] == v2['alt']:
                is_found = True
                fields1[7] = 'PS=B1'
                break

            # Add missing ALT alleles.
            if set(v1['alt']).issubset(v2['alt']) and len(v2['alt']) > len(v1['alt']):
                diff = len(v2['alt']) - len(v1['alt'])
                for allele in v2['alt']:
                    if allele not in v1['alt']:
                        v1['alt'].append(allele)
                if vcf_ad:
                    fields1[9:] = [x + ',0' * diff for x in fields1[9:]]

            # Fix the order of ALT alleles.
            if set(v1['alt']) == set(v2['alt']):
                is_found = True
                mapping = {0: 0}
                for i in range(len(v1['alt'])):
                    mapping[i + 1] = v2['alt'].index(v1['alt'][i]) + 1
                fields1[4] = ','.join(v2['alt'])
                fields1[7] = 'PS=B2'

                # This function returns the fixed GT and AD fields.
                def f(x):
                    gt = x.split(':')[0].split('/')
                    for i in [0, 1]:
                        if gt[i] != '0' and gt[i] != '.':
                            gt[i] = str(mapping[int(gt[i])])
                    if not vcf_ad:
                        return '/'.join(gt)
                    ad1 = x.split(':')[1].split(',')
                    ad2 = [0 for y in ad1]
                    for i in range(len(ad2)):
                        ad2[mapping[i]] = ad1[i]
                    return '/'.join(gt) + ':' + ','.join(ad2)

                fields1[9:] = [f(x) for x in fields1[9:]]
                break

        if not is_found:
            fields1[7] = 'PS=C1'

        conformed_vcf.data.append(fields1)

    return conformed_vcf

def phase_vcf(conformed_vcf, hap_panel, vcf_sep, vcf_ad, out, tr, pd, imp):
    phaseme_vcf = conformed_vcf.copy(["meta", "header", "data"])

    # Temporarily remove the markers that can not be phased statistically.
    for i in reversed(range(len(phaseme_vcf.data))):
        fields = phaseme_vcf.data[i]
        v = parse_vcf_fields(fields)
        if 'PS=C' in fields[7]:
            del phaseme_vcf.data[i]

    if not phaseme_vcf.data:
        logger.warning('Input VCF file will not be phased bacuse there are no eligible markers')
        return conformed_vcf

    if len(phaseme_vcf.data) == 1:
        combined_vcf = conformed_vcf.copy(["meta", "header", "data"])

        for i in range(len(combined_vcf.data)):
            fields = combined_vcf.data[i]
            v = parse_vcf_fields(fields)
            if 'PS=C' not in fields[7]:
                logger.warning('Only marker to be phased: %s' % v['pos'])
                fields[7] = 'PS=D2'
                combined_vcf.data[i] = fields[:9] + [x.replace('/', '|') for x in fields[9:]]
                break
        return combined_vcf

    phaseme_vcf.to_file(f'{out}/phaseme.vcf')

    command = [
        'java', '-Xmx2g', '-jar', f'{pd}/beagle.27Apr20.b81.jar',
        f'gt={out}/phaseme.vcf',
        f'chrom={tr}',
        f'ref={hap_panel}',
        f'out={out}/phased',
        f'impute={str(imp).lower()}'
    ]

    subprocess.call(
        command,
        stdout=open(f'{out}/phased.log', 'a'),
        stderr=open(f'{out}/phased.log', 'a'),
    )

    subprocess.call(['gunzip', f'{out}/phased.vcf.gz'])

    phased_vcf = VCFFile(f'{out}/phased.vcf')
    phased_vcf.read()
    phased_vcf.close()

    combined_vcf = conformed_vcf.copy(["meta", "header"])

    for fields1 in conformed_vcf.data:
        v1 = parse_vcf_fields(fields1)

        # Manually phase the loci where all the samples are homozygous.
        if 'PS=C' in fields1[7]:
            fields1[9:] = [x.replace('./.', '0/0') for x in fields1[9:]]
            if all([x.split(':')[0].split('/')[0] == x.split(':')[0].split('/')[1] for x in fields1[9:]]):
                fields1[9:] = [x.replace('/', '|') for x in fields1[9:]]
                fields1[7] = 'PS=D2'
            combined_vcf.data.append(fields1)
            continue

        is_found = False

        for fields2 in phased_vcf.data:
            v2 = parse_vcf_fields(fields2)
            gt = lambda x: fields2[9:][x].split(':')[0]
            ad = lambda x: fields1[9:][x].split(':')[1]
            idx = list(range(len(fields1[9:])))
            if v1['pos'] == v2['pos'] and v1['ref'] == v2['ref'] and v1['alt'] == v2['alt']:
                is_found = True
                fields1[9:] = [gt(i) + ':' + ad(i) for i in idx] if vcf_ad else [gt(i) for i in idx]
                fields1[7] = 'PS=D1'
                break

        if not is_found:
            logger.warning('Marker filtered by Beagle: %s' % v1['pos'])
            fields1[7] = 'PS=E'

        combined_vcf.data.append(fields1)

    return combined_vcf

def annotate_vcf(combined_vcf, vcf_ad, dt, snpdb, gr):

    annotated_vcf = combined_vcf.copy(["meta", "header", "data"])

    undetected_revertants = [x for x in snpdb if x.rv == 'revertant']

    for fields in annotated_vcf.data:
        v = parse_vcf_fields(fields)

        vi_list = []
        rv_list = []
        so_list = []
        fe_list = []

        for var in v['alt']:
            filtered = [x for x in snpdb if v['pos'] == x.pos and v['ref'] == x.hg and (var == x.var or var == x.wt)]
            if filtered:
                vi = filtered[0].vi
                rv = filtered[0].rv
                so = filtered[0].so
                fe = filtered[0].fe
                undetected_revertants = [x for x in undetected_revertants if x != filtered[0]]

            else:
                vi = 'unknown'
                rv = 'unknown'
                so = 'unknown'
                fe = 'unknown'

            vi_list.append(vi)
            rv_list.append(rv)
            so_list.append(so)
            fe_list.append(fe)

        v['info'].append('VI=' + ','.join(vi_list))
        v['info'].append('RV=' + ','.join(rv_list))
        v['info'].append('SO=' + ','.join(so_list))
        v['info'].append('FE=' + ','.join(fe_list))

        fields[7] = ';'.join(v['info'])

    if dt == "chip":
        pass

    else:
        # Manually add variants that are part of the genome assembly.
        if vcf_ad:
            dat = ["0|0:0,0" for x in annotated_vcf.header[9:]]
            fmt = "GT:AD"

        else:
            dat = ["0|0" for x in annotated_vcf.header[9:]]
            fmt = "GT"

        for snp in undetected_revertants:
            inf = ';'.join([
                'PS=D2',
                f'VI={snp.vi}',
                'RV=revertant',
                f'SO={snp.so}',
                f'FE={snp.fe}'
            ])

            fields = [
                gr['chr'].replace("chr", ""),
                snp.pos,
                snp.id,
                snp.hg,
                snp.wt,
                '.',
                '.',
                inf,
                fmt
            ] + dat

            annotated_vcf.data.append(fields)

        annotated_vcf.data.sort(key = lambda x: int(x[1]))

    return annotated_vcf

def account_vcf(annotated_vcf, vcf_ad):
    accounted_vcf = annotated_vcf.copy(["meta", "header", "data"])

    for fields in accounted_vcf.data:
        v = parse_vcf_fields(fields)

        rv_list = [x for x in v['info'] if 'RV=' in x][0].replace('RV=', '').split(',')
        if 'revertant' not in rv_list:
            continue
        i = rv_list.index('revertant')
        fields[3] = v['alt'][i]
        fields[4] = ','.join([v['ref'] if x == v['alt'][i] else x for x in v['alt']])

        def f(x):
            gt = x.split(":")[0]

            if "|" in gt:
                sep = "|"
            elif "/" in gt:
                sep = "/"
            else:
                raise ValueError("Genotype separator not found")

            field = f"{sep}".join([str(i + 1) if y == '0' else '0' if y == str(i + 1) else y for y in gt.split(sep)])
            if vcf_ad:
                ad = x.split(':')[1].split(',')
                field += ':' + ','.join([ad[i + 1] if y == 0 else ad[0] if y == i + 1 else ad[y] for y in range(len(ad))])
            return field

        fields[9:] = [f(x) for x in fields[9:]]

    return accounted_vcf

def apply_pbe_algorithm(x, v, sample, stardb):
    '''
    This function performs the phase-by-extension (PBE) algorithm. PBE is
    applied to one marker at a time within each sample. PBE will compute four
    scores in the following data strucutre [[A, B], [C, D]]:

        A = 1st allele of GT field (e.g. 0 in 0/1) phased to left (1|0).
        B = 1st allele of GT field (e.g. 0 in 0/1) phased to right (0|1).
        C = 2nd allele of GT field (e.g. 1 in 0/1) phased to left (1|0).
        D = 2nd allele of GT field (e.g. 1 in 0/1) phased to right (0|1).

    Allele with the highest score determines final phasing. These scores
    essentially represent the number of "matching" SNPs for the star allele
    that is defined with the variant of interest, assuming the variant is
    phased to left or right. Here, matching SNPs mean 1) they are used to
    define the specific star allele, 2) they are detected from the sample, and
    3) they were statistically phased to the same haplotype from the previous
    step. If there are more than one star allele defined with the variant of
    interest, PBE will test against all such star alleles and use the highest
    score for the allele.

    Examples:

        GT = 0/1 & Scores = [[0, 0], [0, 0]] -> Result = 0|1 (not "flipped")
        GT = 0/1 & Scores = [[0, 0], [3, 0]] -> Result = 1|0 ("flipped")
        GT = 1/2 & Scores = [[1, 3], [3, 1]] -> Result = 2|1 ("flipped")
    '''

    gt = x.split(':')[0].split('/')

    # If the sample is homozygous, return phased genotype right away.
    if gt[0] == gt[1]:
        return x.replace('/', '|')

    scores = [[0, 0], [0, 0]]

    for i in [0, 1]:

        # If the allele is WT, the score will be 0 automatically.
        if gt[i] == '0':
            continue

        snp1 = SNPAllele()
        snp1.pos = v['pos']
        snp1.wt = v['ref']
        snp1.var = v['alt'][int(gt[i]) - 1]

        stars = [v for k, v in stardb.items() if snp1 in v.core + v.tag]

        for j in [0, 1]:
            for star in stars:
                score = 0
                for snp2 in sample.hap[j].obs:
                    if snp2 in star.core + star.tag:
                        score += 1
                
                # Update with the bigger score.
                if score > scores[i][j]:
                    scores[i][j] = score

    a = scores[0][0]
    b = scores[0][1]
    c = scores[1][0]
    d = scores[1][1]

    if max([a, b]) == max([c, d]):
        if a < b and c > d:
            flip = True
        elif a == b and c > d:
            flip = True
        elif a < b and c == d:
            flip = True
        else:
            flip = False
    else:
        if max([a, b]) > max([c, d]):
            if a > b:
                flip = False
            else:
                flip = True
        else:
            if c > d:
                flip = True
            else:
                flip = False

    if flip:
        result = f'{gt[1]}|{gt[0]}'
    else:
        result = f'{gt[0]}|{gt[1]}'

    if 'AD' in v['format']:
        result = result + ':' + x.split(':')[1]

    he = ','.join([str(x) for x in scores[0] + scores[1]])

    return f'{result}:{he}'

def extend_vcf(accounted_vcf, stardb):
    finalized_vcf = accounted_vcf.copy(["meta", "header", "data"])
    pseudo_biosamples = vcf2biosamples(finalized_vcf, filter=True)
    for fields in finalized_vcf.data:
        v = parse_vcf_fields(fields)
        if 'PS=D' in v['info'][0]:
            continue
        new_data = []
        for i in range(9, len(fields)):
            result = apply_pbe_algorithm(fields[i], v, pseudo_biosamples[i-9], stardb)
            new_data.append(result)
        if not all(new_data):
            logger.warning(f'Marker could not be phased by PBE algorithm: ' % v['pos'])
            continue
        fields[7] = ';'.join(v['info'])
        fields[8] += ':HE'
        fields[9:] = new_data
    return finalized_vcf

def find_candidate_stars(biosamples, gb, stardb, gr):
    for biosample in biosamples:
        biosample.hap[0].start = biosample.hap[1].start = int(gr[f"{gb}_start"])
        biosample.hap[0].end = biosample.hap[1].end = int(gr[f"{gb}_end"])
        f = lambda x: sorted([v for k, v in stardb.items() if set(v.core).issubset(x) and not (v.sv and v.sv not in biosample.sv)], key = lambda x: x.rank)
        hap1_snp = [x for x in biosample.hap[0].obs if x.wt != x.var]
        hap2_snp = [x for x in biosample.hap[1].obs if x.wt != x.var]
        biosample.hap[0].cand = f(hap1_snp)
        biosample.hap[1].cand = f(hap2_snp)
        biosample.dip_cand = f(list(set(hap1_snp + hap2_snp)))

def make_sv_calls(biosamples, tg, out, pd, cg, gdf, dt, cr, gb, sl):
    subprocess.call(["mkdir", f"{out}/ssr"])
    if "chr" in cg or ":" in cg:
        control_gene = "."
    else:
        control_gene = cg
    
    sample_list = ",".join(sl) if sl else "."
    exit_code = subprocess.call([
        "Rscript",
        f"{pd}/sv.R",
        pd,
        out,
        dt,
        tg,
        control_gene,
        cr,
        gdf,
        sample_list,
        gb,
    ])
    if exit_code != 0:
        raise TypeError("Something bad happended during SV detection!")

    for biosample in biosamples:
        with open(f"{out}/ssr/{biosample.name}.txt") as f:
            header = next(f).strip().split("\t")
            fields = next(f).strip().split("\t")
            biosample.sv = [fields[header.index("seq1")], fields[header.index("seq2")]]
            biosample.ssr = fields[header.index("ssr")]

def make_genotype_calls(biosamples, tg, stardb):
    deletion = [v.name for k, v in stardb.items() if v.sv == 'cnv0'][0]
    for biosample in biosamples:
        if biosample.sv == ['no_sv', 'no_sv']:
            biosample.gt = True
        call_sv1(biosample, 'cnv0', deletion, stardb)
        call_cnv3(biosample)
        call_cnv_plus(biosample)
        callers = {
            'cyp2a6': cyp2a6,
            'cyp2b6': cyp2b6,
            'cyp2d6': cyp2d6,
            'cyp2e1': cyp2e1,
            'gstm1': gstm1,
            'gstt1': gstt1,
            'slc22a2': slc22a2,
            'slco1b1': slco1b1,
            'ugt1a4': ugt1a4,
            'ugt2b15': ugt2b15,
            'ugt2b17': ugt2b17
        }
        if tg in callers:
            callers[tg](biosample, stardb)

def plot_profiles(biosamples, out, pd, dp):
    subprocess.call(['mkdir', f'{out}/af'])
    for biosample in biosamples:
        with open(f'{out}/af/{biosample.name}.txt', 'w') as f:
            header = ['pos', 'hap1_al', 'hap2_al', 'hap1_af', 'hap2_af']
            f.write('\t'.join(header) + '\n')
            filtered = [x for x in biosample.hap[0].obs if x.td > 10 and x.het]
            if not filtered:
                continue
            for snp in filtered:
                hap1_snp = [x for x in biosample.hap[0].obs if x.pos == snp.pos][0]
                hap2_snp = [x for x in biosample.hap[1].obs if x.pos == snp.pos][0]
                fields = [snp.pos, hap1_snp.var, hap2_snp.var, str(hap1_snp.af), str(hap2_snp.af)]
                f.write('\t'.join(fields) + '\n')

    subprocess.call(['mkdir', f'{out}/plot'])
    exit_code = subprocess.call(['Rscript', f'{pd}/plot.R', out, str(dp)])
    if exit_code != 0:
        raise TypeError('Something bad happended while plotting profiles!')

def genotype(
        dt: str,
        gb: str,
        tg: str,
        vcf: str,
        out: str,
        cg: Optional[str] = None,
        gdf: Optional[str] = None,
        ref: Optional[str] = None,
        sl: Optional[List[str]] = None,
        dp: bool = False,
        imp: bool = False,
    ) -> None:
    """
    Call star alleles in the target gene from genomic data.

    Args:
        dt (str): Data type (wgs, ts, chip).
        gb (str): Genome build (hg19, hg38).
        tg (str): Target gene.
        vcf (str): VCF file.
        out (str): Output project directory.
        cg (str, optional): Control gene or region.
        gdf (str, optional): GDF file.
        ref (str, optional): Reference VCF file.
        sl (list[str], optional): Sample list.
        dp (bool): Output more detailed plots.
        imp (bool): Impute ungenotyped markers.
    """

    pd = os.path.dirname(os.path.realpath(__file__))
    gene_table = read_gene_table(f"{pd}/gene_table.txt")
    snp_table = read_snp_table(f"{pd}/snp_table.txt", gene_table)
    snpdb = build_snpdb(tg, gb, snp_table)
    star_table = read_star_table(f"{pd}/star_table.txt")
    stardb = build_stardb(tg, gb, star_table, snpdb)
    gr = gene_table[tg]
    tr = gene_table[tg][f"{gb}_region"].replace("chr", "")

    # Log the target region.
    logger.info(f"Target region: chr{tr}")

    # Read the VCF file.
    input_vcf = VCFFile(vcf)
    input_vcf.read(tr, tidy=True)
    input_vcf.close()

    # Log summary data for the VCF file.
    logger.info(f"Number of samples: {len(input_vcf.header) - 9}")
    logger.info(f"Number of markers in target region: {len(input_vcf.data)}")

    if gdf:

        # Determine the control region.
        if not cg:
            raise ValueError("Could not find control gene or region")

        elif "chr" in cg or ":" in cg:
            cr = cg.replace("chr", "")

        else:
            cr = gene_table[cg][f"{gb}_region"].replace("chr", "")

        # Log the control region.
        logger.info(f"Control region: chr{cr}")

        # Assess the GDF file.
        with open(gdf) as f:
            g = [x.replace("Depth_for_", "") for x 
                in f.readline().strip().split("\t")[3:]]

        v = input_vcf.header[9:]

        if len(g) != len(v):
            message = "VCF and GDF files have different sample sizes"
            raise TypeError(message)

        if len(set(g + v)) != len(v):
            message = "VCF and GDF files have different sets of samples"
            raise TypeError(message)

        for i in range(len(v)):
            if v[i] != g[i]:
                message = ("Samples are ordered differently "
                    "between VCF and GDF files")
                raise TypeError(message)

        if dt == "ts" and len(g) < 5:
            message = ("Copy number analysis with TS data "
                "requires at least five samples")
            raise TypeError(message)

    else:
        pass

    # Determine the reference VCF file.
    if ref:
        hap_panel = ref
    else:
        hap_panel = f"{pd}/1kgp/{gb}/{tg}.vcf.gz"

    # Read the reference VCF file.
    ref_vcf = VCFFile(hap_panel)
    ref_vcf.read(tidy=True)
    ref_vcf.close()

    # Make the output project directory.
    os.mkdir(out)

    # Determine the genotype mode.
    vcf_sep, vcf_ad = assess_vcf(input_vcf)

    # Log information about the genotype mode.
    logger.info(f"GT field separator: {vcf_sep}")
    logger.info(f"AD field will be used: {vcf_ad}")

    # Conform genotype fields, filter markers, and add meta data.
    processed_vcf = process_vcf(input_vcf, vcf_sep, vcf_ad, tg, gb)

    # Handle loci with multiple indels that also define star alleles.
    adjusted_vcf = adjust_vcf(processed_vcf, stardb)

    # Conform the input VCF file to match the reference VCF file.
    # Then, statistically phase the input VCF file using the Beagle program.
    if vcf_sep == "|" or not adjusted_vcf.data:
        combined_vcf = conformed_vcf = adjusted_vcf
    else:
        conformed_vcf = conform_vcf(adjusted_vcf, ref_vcf, vcf_ad)
        combined_vcf = phase_vcf(conformed_vcf, hap_panel, vcf_sep, vcf_ad, out, tr, pd, imp)

    # Annotate the input VCF file for various properties.
    annotated_vcf = annotate_vcf(combined_vcf, vcf_ad, dt, snpdb, gr)
    
    # Change the reference sequence from hg19 (or hg38) to the WT allele.
    accounted_vcf = account_vcf(annotated_vcf, vcf_ad)
    
    # Phase the input VCF file again by the haplotype extension algorithm.
    finalized_vcf = extend_vcf(accounted_vcf, stardb)
    finalized_vcf.to_file(f"{out}/finalized.vcf")

    biosamples = vcf2biosamples(finalized_vcf, filter=False)

    if gdf:
        make_sv_calls(biosamples, tg, out, pd, cg, gdf, dt, cr, gb, sl)

    find_candidate_stars(biosamples, gb, stardb, gr)
    make_genotype_calls(biosamples, tg, stardb)
    remove_extra_s1(biosamples)

    # Sort the BioHaplotype objects within each sample.
    for biosample in biosamples:

        if not biosample.gt:
            continue

        a = biosample.hap[0].cand[0].name
        b = biosample.hap[1].cand[0].name
        should_flip = sort_star_names([a, b])[0] == b

        if should_flip:
            biosample.hap[0], biosample.hap[1] = (
                biosample.hap[1], biosample.hap[0])

    predict_phenotypes(biosamples, tg, pd)
    write_result_file(biosamples, out)

    if gdf:
        plot_profiles(biosamples, out, pd, dp)
