# sgRNAble

Brief Outline

## Installation

### Prerequisites

What things you need to install the software and how to install them

* Python3
* Environment Manager (Anaconda is used)


### Installation Guide

Prior to installation,it is best practise to create a new enviroment to store the program and dependencies locally. This setup will create an conda environment with the name sgRNAble and install all required dependencies. 

```
conda create --name sgRNAble python=3.7
conda activate sgRNAble
pip install sgRNAble
conda deactivate
```

End with an example of getting some data out of the system or using it for a little demo

## Quick Run Guide

Ensure that you have a file containing the gene of interest (Target Sequence), the genome of the organism (Genome), and
any additional DNA present. The gene of interest must be present in the genome or the other additional DNA added to the script.

To run the script, ensure you have the following libraries installed:
  - Biopython
  - Tkinter

navigate to the folder in terminal and type in

```
python optimal_guide_finder/guide_finder.py -t data/Fasta_Files/GFP.fasta -g data/Fasta_Files/E_coli_MG1655_genome.fasta data/Fasta_Files/GFP.fasta
```

## Authors
Siddarth Raghuvanshi
Ahmed Abdelmoneim
Avery Noonan
