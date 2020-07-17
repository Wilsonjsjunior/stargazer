from setuptools import setup, find_packages

exec(open("stargazer/version.py").read())

setup(
    name="stargazer",
    version=__version__,
    author='Seung-been "Steven" Lee',
    author_email="sbstevenlee@gmail.com",
    description="A bioinformatics tool for calling star alleles",
    long_description=open("README.rst").read(),
    long_description_content_type="text/x-rst",
    packages=find_packages(),
    entry_points={"console_scripts": ["stargazer=stargazer.__main__:main"]}
)