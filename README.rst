Stargazer
*********

Stargazer is a bioinformatics tool for calling star alleles (haplotypes) 
in pharmacogenes (PGx genes) using data from next-generation 
sequencing (NGS) or single nucleotide polymorphism (SNP) array.

Stargazer is free for academic use, and commercial businesses wanting to 
acquire Stargazer are required to contact license@uw.edu and to purchase a 
commercial license agreement. To download Stargazer, please register 
first by following this 
`link <https://stargazer.gs.washington.edu/stargazerweb/res/form.html>`_.

Installation
============

After the download is complete, move to the Stargazer directory 
and then enter::

    $ python setup.py install

Check if Stargazer is successfully installed::

    $ stargazer -h
    usage: stargazer [-h] [--cg STR] [--gdf FILE] [--ref FILE]
                     [--sl [STR [STR ...]]] [--dp] [--imp]
                     dt gb tg vcf out

    positional arguments:
      dt                    data type (wgs, ts, chip)
      gb                    genome build (hg19, hg38)
      tg                    target gene
      vcf                   VCF file
      out                   output project directory

    optional arguments:
      -h, --help            show this help message and exit
      --cg STR              control gene or region
      --gdf FILE            GDF file
      --ref FILE            reference VCF file
      --sl [STR [STR ...]]  sample list
      --dp                  output more detailed plots
      --imp                 impute ungenotyped markers

Running
=======

Move to the example directory. 
Below uses VDR as the control gene for copy number (CN) analysis::

    $ stargazer \
      wgs \
      hg19 \
      cyp2d6 \
      getrm-cyp2d6-vdr.joint.filtered.vcf \
      ./ex1-getrm-cyp2d6-vdr \
      --gdf getrm-cyp2d6-vdr.gdf \
      --cg vdr

You can provide a custom region as control::

$ stargazer \
  wgs \
  hg19 \
  cyp2d6 \
  getrm-cyp2d6-vdr.joint.filtered.vcf \
  ./ex2-getrm-cyp2d6-custom \
  --gdf getrm-cyp2d6-vdr.gdf \
  --cg chr22:42546883-42551883

Unlike whole genome sequencing data (WGS), target sequencing (TS) data 
require inter-sample normalization for CN analysis. Below example uses 
the population mean during inter-sample normalization::

$ stargazer \
  ts \
  hg19 \
  cyp2d6 \
  hapmap-cyp2d6-vdr.joint.filtered.vcf \
  ./ex3-hapmap-cyp2d6-vdr \
  --gdf hapmap-cyp2d6-vdr.gdf \
  --cg vdr

You may indicate known reference samples without structural variation.
Below uses the mean of indicated samples instead of the population mean::

$ stargazer \
  ts \
  hg19 \
  cyp2d6 \
  hapmap-cyp2d6-vdr.joint.filtered.vcf \
  ./ex4-hapmap-cyp2d6-vdr-list \
  --gdf hapmap-cyp2d6-vdr.gdf \
  --cg vdr \
  --sl 133419 133420 133421 133423 133425

Below runs in VCF only mode for hg19 data::

$ stargazer \
  wgs \
  hg19 \
  cyp3a5 \
  getrm-cyp3a5-hg19.joint.filtered.vcf \
  ex5-getrm-cyp3a5-vcfonly-hg19

Run with hg38 data::

$ stargazer \
  wgs \
  hg38 \
  cyp3a5 \
  getrm-cyp3a5-hg38.joint.filtered.vcf \
  ex6-getrm-cyp3a5-vcfonly-hg38

Run with chip data::

$ stargazer \
  chip \
  hg19 \
  cyp3a5 \
  rok-cyp3a5.vcf \
  ex7-rok-cyp3a5

Run with imputation of ungenotyped markers::

$ stargazer \
  chip \
  hg19 \
  cyp3a5 \
  rok-cyp3a5.vcf \
  ex8-rok-cyp3a5 \
  --imp