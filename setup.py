import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README_PYPI.md").read_text()

setup(
    name='epilogos',
    version='1.0.0',
    description="Information-theoretic navigation of multi-tissue functional genomic annotations",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/meuleman/epilogos",
    authors=['Wouter Meuleman', 'Jacob Quon', 'Alex Reynolds', 'Eric Rynes'],
    license='GNU General Public License v3.0',
    packages=['epilogos'],
    include_package_data=True,
    install_requires=['cython', 'click', 'numpy', 'pandas', 'pyranges', 'scipy', 'matplotlib', 'statsmodels'],
    entry_points={
        'console_scripts': [
            'epilogos = epilogos:main',
        ],
    }
)
