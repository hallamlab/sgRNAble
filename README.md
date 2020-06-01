# sgRNAble

Abstract

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

## Quick Run Guide

Ensure that you have a file containing the gene of interest (Target Sequence), the genome of the organism (Genome), and
any additional DNA present. The gene of interest must be present in the genome or the other additional DNA added to the script.

navigate to the folder in terminal and type in

```
pip install .

sgrnable -t data/Fasta_Files/GFP.fasta -g data/Fasta_Files/E_coli_MG1655_genome.fasta data/Fasta_Files/GFP.fasta
```

## Authors
* [Siddarth Raghuvanshi](https://github.com/Siddarth-Raghuvanshi)
* [Ahmed Abdelmoneim](https://github.com/AhmedAbdelmoneim)
* Avery Noonan

## Contact

Need something? Send me an email at Raghuvanshi.Siddarth@gmail.com
