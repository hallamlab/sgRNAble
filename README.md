# CRISPR_Guide_RNA

Uses the model developed by Farasat/Salis to run a script that finds the best guide RNAs
for a given gene (Target Sequence), Genome and any additional DNA sequences. The model is based on the
information based on CAS9 and would be different for other CRISPR systems

## Getting Started

These instructions will get you a copy of the project up and running on your local machine

### Prerequisites

What things you need to install the software and how to install them

```
python3
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
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
python optimal_guide_calc/optimal_guide.py -t data/Fasta_Files/GFP.fasta -g data/Fasta_Files/E_coli_MG1655_genome.fasta data/Fasta_Files/GFP.fasta
```

Follow the prompts on screen.

## Running the tests

Explain how to run the automated tests for this system

### And coding style tests

Ensure Pylint is installed:

You can either set-up Pylint in your development enviornment directly or run through the command line

```
# in the project's root directory run:

pylint --max-line-length=120 optimal_guide_calc
```

Make sure pylint passes with no errors or warning before commiting your code.

In addition make sure any added code is documented, function docstrings should follow the following format:

```
    """
    [summary]

    Arguments:
        arg1 {[type]} -- [description]
        arg2 {[type]} -- [description]

    Raises:
        exceptionType: [description]

    Returns:
        [type] -- [description]
    """
```

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.



## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
