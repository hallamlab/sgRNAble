import numpy as np
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
from Cas9_Calculator import *
from argparse import RawTextHelpFormatter
import time
from Azimuth_Model import model_comparison

# Returns the upper case sequences as strings from the files given as arguments. Also combines the various genome sequences
def Get_Sequence(args):

    #Reads the file using biopython and creates a object called target
    Target_Dict = SeqIO.to_dict(SeqIO.parse(args.target_sequence, args.target_sequence.split('.')[1]))

    for name in Target_Dict:
        Target_Dict[name] = Target_Dict[name].seq.upper()

    #Reads the Genome files using biopython and combines them into one genome object
    Genome = SeqIO.read(args.genome_sequence[0], args.genome_sequence[0].split('.')[1])
    for i in range(1,len(args.genome_sequence)):
        Genome  = Genome + SeqIO.read(args.genome_sequence[i], args.genome_sequence[i].split('.')[1])

    return Target_Dict, Genome.seq.upper()

#Find the Guide RNAs in a Sequence
def PAM_Finder(Sequence,args):

    #initalize arguements into variables for the user
    PAM = "GG"
    guide_RNA_length = args.guide_length
    Azimuth_Distance = 24 # from the from of the NGG sequence
    cut_list = args.cut
    if cut_list == None:
        cut_list = []

    Locations = []

    Sequence = str(Sequence)
    Position = 0
    while Position < len(Sequence):
        i = Sequence[Position:].find(PAM)
        #Finds the location of the next cut argument which might be an issue
        potential_guide_location = Position + i - guide_RNA_length
        discard = any(cutsite in Sequence[potential_guide_location: potential_guide_location + guide_RNA_length] for cut_site in cut_list)
        if i < 0:
            break

        if potential_guide_location > (Azimuth_Distance - args.guide_length) and not discard:
            Locations.append(potential_guide_location - 1)

        Position = Position + i + 1

    return Locations

def Guide_Selection(Target_Dict, args):

    #Obtain the Guide RNAs from the Target Sequence
    guide_list = {}

    #Default purpose, the tool looks for all avalible guides possibly found in the sequence given by the user
    if args.purpose == "d":
        for gene in Target_Dict:
            #Get guides on the positive strand
            guide_list[gene] = []
            Locations = PAM_Finder(Target_Dict[gene], args)
            strand_array = ["Positive"] * len(Locations)

            #Get guides on the negative strand
            NegLocations = PAM_Finder(Seq.reverse_complement(Target_Dict[gene]), args)
            Sequence_length = len(Target_Dict[gene])
            NegLocations = [Sequence_length - x for x in NegLocations]
            negative_array = ["Negative"] * len(NegLocations)

            #Combine the information into a single array for each gene
            guides = [str(Target_Dict[gene][loc-4:loc+args.guide_length+6]) for loc in Locations]
            guides.extend([str(Seq.reverse_complement(Target_Dict[gene][loc-6-args.guide_length:loc+4])) for loc in NegLocations ])
            Locations.extend(NegLocations)
            strand_array.extend(negative_array)

            #Run the model through the Azimuth Model
            predictions = model_comparison.predict(np.array(guides))
            guides = [x for y,x in sorted(zip(predictions,guides), reverse=True)][:args.azimuth_cutoff]
            Locations =  [x for y,x in sorted(zip(predictions,Locations), reverse=True)][:args.azimuth_cutoff]
            strand_array =  [x for y,x in sorted(zip(predictions,strand_array), reverse=True)][:args.azimuth_cutoff]
            guides = [guide[4:24] for guide in guides]
            guide_list[gene].extend([guides, Locations, strand_array])

    #only runs the negative strand as CRISPRi works better on negative strands
    elif args.purpose == "i":
        for gene in Target_Dict:
            Locations = PAM_Finder(Seq.reverse_complement(Target_Dict[gene]), args)
            Sequence_length = len(Target_Dict[gene])
            NegLocations = [Sequence_length - x for x in Locations]
            guide_list[gene].append(NegLocations)
    elif args.purpose == "s":
        #Run through Prodigal

        #am only running the negative side as prodigal gives you the positive side of the gene and looks for the gene on both strands (I'm sure of this, but I should verify it)
        for gene in Target_Dict:
            Locations = PAM_Finder(Seq.reverse_complement(Target_Dict[gene]), args)
            Sequence_length = len(Target_Dict[gene])
            NegLocations = [Sequence_length - x for x in Locations]
            guide_list[gene].append(NegLocations)
    elif args.purpose == "g":
        pass
    else:
        #cull all guides which have a location larger than the user given cutoff as this is the activation case
        for gene in Target_Dict:
            guide_list[gene] = []
            PosLocations = PAM_Finder(Target_Dict[gene], args)
            PosLocations = [x for x in PosLocations if x < int(args.purpose[:-1])]
            guide_list[gene].append(PosLocations)

            Locations = PAM_Finder(Seq.reverse_complement(Target_Dict[gene]), args)
            Sequence_length = len(Target_Dict[gene])
            NegLocations = [Sequence_length - x for x in Locations]
            NegLocations = [x for x in NegLocations if x < int(args.purpose[:-1]) ]
            guide_list[gene].append(NegLocations)

    return guide_list

def main():

    #Parser to get the files listed in the arguments
    parser = argparse.ArgumentParser(description="""This program helps you to find all possible guide RNAs that will \ntarget the gene. Then using the model created by Salis Lab, \nyou can see the off target effects for the each possible guide.""",
                                     formatter_class=RawTextHelpFormatter)

    #Parsers to add arguements.
    parser.add_argument("-t", "--target_sequence", required=True,
                        help= "The Gene Sequence of Interest (Fasta or Genebank)")
    parser.add_argument("-g", "--genome_sequence", required=True, nargs = '+',
                        help= "The Genome of the organism, if targeting a plasmid, make sure to \n include it as well (Fasta or Genebank)")
    parser.add_argument("-a", "--azimuth_cutoff", required=False, default = 10,
                        help= "How many guides should pass from azimuth screening, the guides are passed based on descending azimuth prediction score")
    parser.add_argument("-l", "--guide_length", required = False, default = 20,
                        help = "Length of the guide RNA sequence")
    parser.add_argument("-p", "--purpose", required=False, default = "d",
                        help= " i: CRISPR interference on gene \n ###a: CRISPR activation on gene, enter the number of base pair from start you would want \n s: CRISPRi screening from contigs (genes found via prodigal) \n g: guide binding strength calculator \n Leave blank to see all possible guides and off target effects from your sequence")
    parser.add_argument("-c", "--cut", required=False, nargs = "+",
                        help= "Sequences to avoid in guides (i.e restriction enzyme sites)")

    #Creating a variable to make the values easily accessible
    args = parser.parse_args()

    #Get the sequences in a Seq format from user fasta or genebank files
    Target_Dict, Genome = Get_Sequence(args)

    ref_record = SeqRecord(Genome, id="refgenome", name ="reference", description ="a reference background")
    ref_record = ref_record + ref_record.reverse_complement()
    SeqIO.write(ref_record, "Run_Genome_Plus_RC", "fasta")

    #Select the guides based on the purpose and the azimuth model
    guide_list = Guide_Selection(Target_Dict, args)

    #Build the model
    __start = time.time()
    Cas9Calculator=clCas9Calculator('Run_Genome_Plus_RC')
    #if args.purpose == "g":
        #different target guides
    sgRNA_Created = sgRNA(guide_list, Cas9Calculator)
    __elasped = (time.time() - __start)
    print("Time Model Building: {:.2f}".format(__elasped))

    #Run the model
    __start = time.time()
    sgRNA_Created.run()
    __elasped = (time.time() - __start)
    print("Time model calculation: {:.2f}".format(__elasped))

if __name__ == "__main__":
    main()
