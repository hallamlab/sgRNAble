import numpy as np
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import Seq
from Cas9_Calculator import *
from File_GUI import *


def Get_Sequence():

    Root = Tk()
    Program = GUI(Root)
    Program.run()
    Root.mainloop()
    i = 0
    Genome_Created = False


    for i in range(len(Program.Filetypes)):
        if not(Program.Filetypes[i] == 0):
            if(Program.Files[1][i] == "TARGET"):
                Target = SeqIO.read(Program.Files[0][i], Program.Filetypes[0].get().lower())
            elif not(Genome_Created):
                Genome = SeqIO.read(Program.Files[0][i], Program.Filetypes[i].get().lower())
            else:
                Genome  = Genome + SeqIO.read(Program.Files[0][i], Program.Filetypes[1].get().lower())
        else:
            pass

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
  j = 0 #Variable for limiting time spent searching Genome
  while True:
    i = Temp_Sequence.find(PAM)
    if(i == -1):
        break
    if(j > 10000):
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
    j = j+1

  return Guide_RNAs,Location,Strand

#Combine the Coding and Template Strands into a single strand
def CombinetoStr (Template_Guides, Coding_Guides):
  Guides = []

  for i in range (len(Template_Guides)):
    if (i < len(Template_Guides)):
      Guides.append(str(Template_Guides[i]))

  for i in range (len(Coding_Guides)):
    if (i < len(Coding_Guides)):
      Guides.append(ReverseComplement(str(Coding_Guides[i])))

  return Guides

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
