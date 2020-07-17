from setuptools import setup, find_packages

exec(open("stargazer/version.py").read())

setup(
    name="stargazer",
    version=__version__,
    author='Seung-been "Steven" Lee',
    author_email="sbstevenlee@gmail.com",
    description="A tool for calling star alleles from pharmacogenomic data",
    long_description=open("README.rst").read(),
    long_description_content_type="text/x-rst",
    packages=find_packages(),
    package_data={
        "stargazer": [
            "plot.R",
            "sv.R",
            "sv_table.txt",
            "star_table.txt",
            "snp_table.txt",
            "phenotype_table.txt",
            "gene_table.txt",
            "beagle.27Apr20.b81.jar",
        ],
        "stargazer.1kgp": [
            "hg19/*.vcf.gz",
            "hg38/*.vcf.gz",
        ]
    },
    entry_points={"console_scripts": ["stargazer=stargazer.__main__:main"]}
)