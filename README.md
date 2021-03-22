# sgRNAble

Tool for high-throughput design of sgRNA libraries targeting selected genes or whole genomes, while considering both on-target binding potential and off-target effects of a given sgRNA in a user-defined genome.

## Installation

### Prerequisites

What things you need to install the software and how to install them

* Python3
* Environment Manager (Anaconda is used here)


### Installation Guide

Prior to installation,it is best practise to create a new enviroment to store the program and dependencies locally. This setup will create an conda environment with the name sgRNAble and install all required dependencies. Start this process by navigating to the path of the github download(inside the folder).

```
cd PATH/TO/sgRNAble
conda create --name sgRNAble python=3.7
conda activate sgRNAble
pip install .
conda deactivate
```

In the future, the program can be run by activating the python env and running the program.

```
conda activate sgRNAble
sgrnable -t TARGET_FILE -g GENOME_FILE
conda deactivate
```

### Running Tests
To run package tests, cd into the project directory and run the following command. This installs any missing dependencies and runs package tests.

```
python setup.py test
```

## Quick Run Guide

Ensure that you have a file containing the gene of interest (Target Sequence), the genome of the organism (Genome), and
any additional DNA present. The gene of interest must be present in the genome or the other additional DNA added to the script.

To start a test run targeting GFP in E.Coli genome, navigate to repository root and run the following command:

```
pip install .

sgrnable -t tests/data/gfp.fasta -g tests/data/ecoli_genome.fasta tests/data/gfp.fasta -th 4
```

## Distribution Guide

For pushing tool to PyPi and Conda

### Versioning
Update version number in setup.py. Number is in the format major.minor.patch
* for small updates increment patch number
* for minor features increment minor number
* for major features increment major number

### PyPi
1. Build source distribution by running ```python setup.py sdist```
2. Ensure twine is installed ```pip install twine```
3. Push source to Pypi by running ```twine upload dist/*```

## Authors
* [Siddarth Raghuvanshi](https://github.com/Siddarth-Raghuvanshi)
* [Ahmed Abdelmoneim](https://github.com/AhmedAbdelmoneim)
* [Avery Noonan](https://github.com/Noonanav)

## Contact

Need something? Send me an email at Raghuvanshi.Siddarth@gmail.com

## References

Farasat, I., & Salis, H. M. (2016). A Biophysical Model of CRISPR/Cas9 Activity for Rational Design of Genome Editing and Gene          Regulation. _PLOS Computational Biology_, 12(1), e1004724. doi:10.1371/journal.pcbi.1004724

Doench, J. G., Fusi, N., Sullender, M., Hegde, M., Vaimberg, E. W., Donovan, K. F., . . . Root, D. E. (2016). Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. _Nature biotechnology_, 34, 184. doi:10.1038/nbt.3437
