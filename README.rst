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