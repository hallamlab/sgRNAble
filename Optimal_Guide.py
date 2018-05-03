import numpy
from Bio import SeqIO
from tkinter import Tk
from tkinter import messagebox
from tkinter.filedialog import askopenfilename


messagebox.showinfo("Target Sequence", "Please select a Target Sequence Fasta File")
Tk().withdraw()
Target_Seq_Filename  = askopenfilename()
print(Target_Seq_Filename)

Genome_Filename = "Fasta_Files/E_coli_MG1655_Genome.fasta"
Target = SeqIO.read(Target_Seq_Filename, "fasta")
Genome = SeqIO.read(Genome_Filename, "fasta")

#Convert Sequences to uppercase
Target_Seq = Target.seq.upper()
Genome_Seq = Genome.seq.upper()

#Length of the Guide RNA desired
Guide_RNA_length = 20

#Find the Guide RNAs in a Sequence
def PAM_Finder(Sequence, PAM, strand):
  Guide_RNAs = []

  if (strand == 1):
    Direction = 1
  else:
    Direction = -1
  Position = 0
  Temp_Sequence = Sequence
  while True:
    i = Temp_Sequence.find(PAM)
    if(i == -1):
      break
    Position = Position + i + 2
    if(Position > Guide_RNA_length):
      if(Direction > 0):
        Guide_RNAs.append(Sequence[Position-23:Position])
      if(Direction < 0):
        Guide_RNAs.append(Sequence[Position-2:Position+21])
    Temp_Sequence = Temp_Sequence[i+2:]

  return Guide_RNAs

#Combine the Coding and Template Strands into a single strand
def Combine (Template_Guides, Coding_Guides):
  Guides = []

  for i in range (len(Template_Guides)):
    if (i < len(Template_Guides)):
      Guides.append(str(Template_Guides[i]))

  for i in range (len(Coding_Guides)):
    if (i < len(Coding_Guides)):
      Guides.append(reverse_complement(str(Coding_Guides[i])))

  return Guides

def reverse_complement(nucleotide_sequence):
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

#Obtain the Guide RNAs from the Target Sequence
T_Guides_GG = PAM_Finder(Target_Seq, "GG",1)
T_Guides_CC = PAM_Finder(Target_Seq, "CC", 2)

Target_Guides = Combine(T_Guides_GG, T_Guides_CC)
print(Target_Guides)





#Obtain the all possible off target pam sites in the genome
