import pathlib
from setuptools import setup, find_packages
from epilogos.__init__ import __version__

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.txt").read_text()

# Requirements kept in one location
REQUIREMENTS = (HERE / "requirements.txt").read_text()
install_requirements = REQUIREMENTS.splitlines()

setup(
    name="epilogos",
    version=__version__,
    author="Wouter Meuleman, Jacob Quon, Alex Reynolds, Eric Rynes",
    author_email="wouter@meuleman.org",
    description="Information-theoretic navigation of multi-tissue functional genomic annotations",
    long_description=README,
    long_description_content_type="text/x-rst",
    url="https://github.com/meuleman/epilogos",
    license="LICENSE.txt",
    packages=find_packages("."),
    scripts=["bin/preprocess_data_ChromHMM.sh", "bin/download_example_data.sh"],
    include_package_data=True,
    install_requires=install_requirements,
    entry_points={
        "console_scripts": [
            "epilogos = epilogos.run:main",
            "simsearch = epilogos.similaritySearch_run:main",
            "plotregion = epilogos.plotRegion:main"
        ],
    }
)