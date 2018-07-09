import numpy as np
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
from Cas9_Calculator import *
from argparse import RawTextHelpFormatter


#Parser to get the files listed in the arguments

#Creating the parser
parser = argparse.ArgumentParser(description="""This program helps you to find all possible guide RNAs that will \ntarget the gene. Then using the model created by Salis Lab, \nyou can see the off target effects for the each possible guide.""", formatter_class=RawTextHelpFormatter)

#Parsers to add arguements.
parser.add_argument("-t", "--target_sequence", required=True, help= "The Gene Sequence of Interest (Fasta or Genebank)")
parser.add_argument("-g", "--genome_sequence", required=True, nargs = '+', help= "The Genome of the organism, if targeting a plasmid, make sure to \n include it as well (Fasta or Genebank)")
parser.add_argument("-k", "--kind", required=False, help= " a: CRISPR activation \n i: CRISPR interference \n Leave Blank to see all possible off target effects of your sequence")

#Creating a variable to make the values easily accessible
args = parser.parse_args()

# Returns the upper case sequences as strings from the files given as arguments. Also combines the various genome sequences
def Get_Sequence():

    Target = SeqIO.read(args.target_sequence, args.genome_sequence[0].split('.')[1]) #Reads the file using biopython and creates a object called target

    #Reads the Genome files using biopython and combines them into one genome object
    Genome = SeqIO.read(args.genome_sequence[0], args.genome_sequence[0].split('.')[1])
    for i in range(1,len(args.genome_sequence)):
        Genome  = Genome + SeqIO.read(args.genome_sequence[i], args.genome_sequence[0].split('.')[1])

    return Target.seq.upper(), Genome.upper()

#Length of the Guide RNA desired
Guide_RNA_length = 20

#Find the Guide RNAs in a Sequence
def PAM_Finder(Sequence, PAM, Direction):
  Guide_RNAs = []
  Location = []
  Strand = []

  Position = 0
  Temp_Sequence = Sequence
  while True:
    i = Temp_Sequence.find(PAM)
    if(i == -1):
        break
    Position = Position + i + 2
    if(Position > Guide_RNA_length):
        if(Direction > 0):
            Location.append(Position - 2)
            Strand.append(Direction)
            Guide_RNAs.append(Sequence[Position-23:Position-3])
        if(Direction < 0):
            Location.append(Position + 1)
            Strand.append(Direction)
            Guide_RNAs.append(Sequence[Position+1:Position+21])
    Temp_Sequence = Temp_Sequence[i+2:]

  return Guide_RNAs,Location,Strand

#Combine the Coding and Template Strands into a single strand and converts that into a string
def CombinetoStr (Template_Guides, Coding_Guides):
  Guides = []

  for i in range (len(Template_Guides)):
    if (i < len(Template_Guides)):
      Guides.append(str(Template_Guides[i]))

  for i in range (len(Coding_Guides)):
    if (i < len(Coding_Guides)):
      Guides.append(ReverseComplement(str(Coding_Guides[i])))

  return Guides

#Reverse a strand
def ReverseComplement(nucleotide_sequence):
  comp = []
  for c in nucleotide_sequence:
    if c == 'A' or c == 'a':
      comp.append('T')
    if c == 'G' or c == 'g':
      comp.append('C')
    if c == 'U' or c == 'u' or c == 'T' or c == 't':
      comp.append('A')
    if c == 'C' or c == 'c':
      comp.append('G')
  rev_comp = ''.join(reversed(comp))
  return rev_comp

Target_Seq, Genome = Get_Sequence()

Genome = Genome + Genome.reverse_complement()
SeqIO.write(Genome, "Total_Genome_Plus_RC", "fasta")

#Obtain the Guide RNAs from the Target Sequence
T_Guides_Pos, Position_Pos, Direction_Pos = PAM_Finder(Target_Seq, "GG",1)
T_Guides_Neg, Position_Neg, Direction_Neg = PAM_Finder(Target_Seq, "CC", -1)

#Combine the information from the Guide RNAs into one single array.
Position_List = Position_Pos + Position_Neg
Direction_List = Direction_Pos + Direction_Neg
Guide_Info = np.vstack((Position_List, Direction_List)).T

#Combine the two guides
Target_Guides = CombinetoStr(T_Guides_Pos, T_Guides_Neg)

#Send Data to the model
Cas9Calculator=clCas9Calculator(['Total_Genome_Plus_RC'])
sgRNA1 = sgRNA(Target_Guides,Guide_Info, Cas9Calculator)
sgRNA1.run()
