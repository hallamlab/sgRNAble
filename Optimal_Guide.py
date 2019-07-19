import numpy as np
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
from Cas9_Calculator import *
from argparse import RawTextHelpFormatter
import time

# Returns the upper case sequences as strings from the files given as arguments. Also combines the various genome sequences
def Get_Sequence(args):

    #Reads the file using biopython and creates a object called target
    Target_Dict = SeqIO.to_dict(SeqIO.parse(args.target_sequence, args.genome_sequence[0].split('.')[1]))

    for name in Target_Dict:
        Target_Dict[name] = Target_Dict[name].seq.upper()

    #Reads the Genome files using biopython and combines them into one genome object
    Genome = SeqIO.read(args.genome_sequence[0], args.genome_sequence[0].split('.')[1])
    for i in range(1,len(args.genome_sequence)):
        Genome  = Genome + SeqIO.read(args.genome_sequence[i], args.genome_sequence[i].split('.')[1])

    return Target_Dict, Genome.upper()

#Find the Guide RNAs in a Sequence
def PAM_Finder(Sequence, PAM, guide_RNA_length, Discard):

    Locations = []

    Sequence = str(_Sequence)
    Position = 0
    while Position < len(Sequence):
        i = Sequence[Position:].find(PAM)
        j = Sequence[Position:].find(Discard)
        if i < 0:
            break
        # Doesn't add any guides which contain a site from the cut arguement
        if i - j > 22:
            continue
        if(Position + i - 22 > 0)
            Locations.append(Position + i - 22)

        Position = Position + i + len(PAM)

    return Locations

 def main():

    #Parser to get the files listed in the arguments
    parser = argparse.ArgumentParser(description="""This program helps you to find all possible guide RNAs that will \ntarget the gene. Then using the model created by Salis Lab, \nyou can see the off target effects for the each possible guide.""",
                                     formatter_class=RawTextHelpFormatter)

    #Parsers to add arguements.
    parser.add_argument("-t", "--target_sequence", required=True,
                        help= "The Gene Sequence of Interest (Fasta or Genebank)")
    parser.add_argument("-g", "--genome_sequence", required=True, nargs = '+',
                        help= "The Genome of the organism, if targeting a plasmid, make sure to \n include it as well (Fasta or Genebank)")
    parser.add_argument("-p", "--pam_Sequence", required =False, default = "GG",
                        help= "The PAM sequence you want the code to look for \nLeave blank to use NGG")
    parser.add_argument("-l", "--guide_length", required = False, default = 20,
                        help = "Length of the guide RNA sequence")
    parser.add_argument("-a", "--aim", required=False, default = "d",
                        help= " i: CRISPR interference on gene \n ###a: CRISPR activation on gene, enter the number of base pair from start you would want \n s: CRISPRi screening from contigs \n g: Guide Binding Strength \n Leave blank to see all possible guides and off target effects from your sequence")
    parser.add_argument("-c", "--cut", required=False, nargs = "+"
                        help= "Sequences to avoid in guides (i.e restriction enzyme sites)")

    #Creating a variable to make the values easily accessible
    args = parser.parse_args()

    #Length of the Guide RNA desired
    Guide_RNA_length = args.guide_length

    Target_Dict, Genome = Get_Sequence(args)

    ref_record = SeqRecord(Genome, id="refgenome", name ="reference", description ="a reference background")
    SeqIO.write(ref_record, "Run_Genome_Plus_RC", "fasta")

    #Obtain the Guide RNAs from the Target Sequence
    guide_list = {}

    if args.aim == "d":
        for gene in Target_Dict:
            guide_list(gene) = []
            PosLocations = PAM_Finder(Target_Dict[gene], args.pam_Sequence, Guide_RNA_length)
            guide_list[gene].append(PosLocations)

            Locations = PAM_Finder(Seq.reverse_complement(Target_Dict[gene]), args.pam_Sequence, Guide_RNA_length)
            Sequence_length = len(Target_Dict[gene])
            NegLocations = [Sequence_length - x for x in Locations]
            guide_list[gene].append(NegLocations)
    elif args.aim == "i":
        #only runs the negative strand as CRISPRi works wayy better on negative strands
        for gene in Target_Dict:
            Locations = PAM_Finder(Seq.reverse_complement(Target_Dict[gene]), args.pam_Sequence, Guide_RNA_length)
            Sequence_length = len(Target_Dict[gene])
            NegLocations = [Sequence_length - x for x in Locations]
            guide_list[gene].append(NegLocations)
    elif args.aim == "s":
        #Run through Prodigal

        #am only running the negative side as prodigal gives you the positive side of the gene and looks for the gene on both strands
        for gene in Target_Dict:
            Locations = PAM_Finder(Seq.reverse_complement(Target_Dict[gene]), args.pam_Sequence, Guide_RNA_length)
            Sequence_length = len(Target_Dict[gene])
            NegLocations = [Sequence_length - x for x in Locations]
            guide_list[gene].append(NegLocations)
    elif args.aim == "g":
        continue
    else:
        #cull all guides which have a location larger than the user given cutoff as this is the activation case
        for gene in Target_Dict:
            guide_list(gene) = []
            PosLocations = PAM_Finder(Target_Dict[gene], args.pam_Sequence, Guide_RNA_length)
            PosLocations = [x for x in PosLocations if x < int(args.aim[:-1])]
            guide_list[gene].append(PosLocations)

            Locations = PAM_Finder(Seq.reverse_complement(Target_Dict[gene]), args.pam_Sequence, Guide_RNA_length)
            Sequence_length = len(Target_Dict[gene])
            NegLocations = [Sequence_length - x for x in Locations]
            NegLocations = [x for x in NegLocations if x < int(args.aim[:-1]) ]
            guide_list[gene].append(NegLocations)

    if args.aim == "g"
        # Test the effect of running a guide through the SeqIO dictionary thing.
        # Run the model and then just quit()

    #Build the model
    __start = time.time()
    Cas9Calculator=clCas9Calculator(['Total_Genome_Plus_RC'])
    if args.aim == "g":
        #different target guides
    sgRNA_Created = sgRNA(Target_Guides,Guide_Info, Cas9Calculator)
    __elasped = (time.time()- __start)
    print("Time Model Building: {:.2f}".format(__elasped))

    #Run the model
    __start = time.time()
    sgRNA_Created.run()
    __elasped = (time.time() - __start)
    print("Time model calculation: {:.2f}".format(__elapsed))

if __name__ == "__main__":
    main()
