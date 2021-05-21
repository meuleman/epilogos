import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.txt").read_text()

setup(
    name="Epilogos",
    version="fshjfhkjfadsghjkfs",
    authors=["Wouter Meuleman", "Jacob Quon", "Alex Reynolds", "Eric Rynes"],
    description="Information-theoretic navigation of multi-tissue functional genomic annotations",
    long_description=README,
    long_description_content_type="text/x-rst",
    url="https://github.com/meuleman/epilogos",
    license="LICENSE.txt",
    packages=["epilogos"],
    scripts=["bin/preprocess_data_ChromHMM.sh"],
    include_package_data=True,
    install_requires=[
        "cython == 0.29.23",
        "statsmodels == 0.12.0",
        "scipy == 1.5.2",
        "numpy == 1.19.2",
        "matplotlib == 3.3.2",
        "click == 7.1.2",
        "pandas == 1.1.3",
        "pyranges == 0.0.97",
        ],
    entry_points={
        "console_scripts": [
            "epilogos = epilogos.__main__:main",
        ],
    }
)
