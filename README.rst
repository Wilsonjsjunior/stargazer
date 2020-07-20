Stargazer
*********

Stargazer is a bioinformatics tool for calling star alleles (haplotypes) 
in pharmacogenes (PGx genes) using data from next-generation 
sequencing (NGS) or single nucleotide polymorphism (SNP) array. For more 
information, please visit the official Stargazer webpage: 
https://stargazer.gs.washington.edu/stargazerweb/.

Stargazer is free for academic use, and commercial businesses wanting to 
acquire Stargazer are required to contact license@uw.edu and to purchase a 
commercial license agreement. To download Stargazer, please register 
first by following this 
`link <https://stargazer.gs.washington.edu/stargazerweb/res/form.html>`_.

Author
======

Stargazer was developed by Seung-been Lee (he goes by "Steven") during 
his PhD in the Nickerson lab at the University of Washington. Steven 
graduated in June of 2019 and now works in industry, but he's still in 
charge of developing and maintaining Stargazer.

Citation
========

If you use Stargazer in a published analysis, please report the program 
version and cite the appropriate article. The most recent reference for 
Stargazer's genotyping algorithm is:

Lee et al., 2019. Calling star alleles with Stargazer in 28 pharmacogenes 
with whole genome sequences. Clinical Pharmacology & Therapeutics. 
DOI: https://doi.org/10.1002/cpt.1552.

The original genotyping pipeline using Stargazer is described in:

Lee et al., 2018. Stargazer: a software tool for calling star alleles 
from next-generation sequencing data using CYP2D6 as a model. 
Genetics in Medicine. DOI: https://doi.org/10.1038/s41436-018-0054-0.

Installation
============

After the download is complete, move to the Stargazer directory 
and then enter::

    $ python setup.py install

Check if Stargazer is successfully installed by entering::

    $ stargazer -h

To give::

    usage: stargazer [-h] [--version] [--cg STR] [--gdf FILE] [--ref FILE]
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
      --version             print the Stargazer version number and exit
      --cg STR              control gene or region
      --gdf FILE            GDF file
      --ref FILE            reference VCF file
      --sl [STR [STR ...]]  sample list
      --dp                  output more detailed plots
      --imp                 impute ungenotyped markers

Running
=======

Stargazer supports many different genotyping modes. In this section, we 
will show some examples using real data. To do this, first move to the 
example directory.

Example 1
---------

Below will genotype 70 Coriell DNA samples for the CYP2D6 gene using their
whole genome sequencing (WGS) data (taken from Lee et al., 2019). 
The VDR gene will be used as the control gene for copy number (CN) analysis::

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