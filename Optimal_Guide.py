import numpy as np
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import Seq
from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter import messagebox
from Cas9_Calculator import *


def Get_Files (Title, Repeat):

    root = Tk()
    Message = "Please select the " + Title + " file"
    More_Seqs = True

    # Add a grid
    mainframe = Frame(root)
    mainframe.grid(column=0,row=0, sticky=(N,W,E,S) )
    mainframe.columnconfigure(0, weight = 1)
    mainframe.rowconfigure(0, weight = 1)
    mainframe.pack(pady = 100, padx = 100)

    if(Repeat):
        Message = "Please select an " + Title + " file"
        More_Seqs = messagebox.askyesno( Title,"Would you like to add Additional Sequences to the genome \n (Such as Plasmids)")
        if not(More_Seqs):
            return 1

    messagebox.showinfo(Title, "Please select the file type")

    # Create a Tkinter variable
    tkvar = StringVar(root)

    # Dictionary with options
    Options = [ "Fasta", "Genbank",]
    tkvar.set('Fasta') # set the default option

    popupMenu = OptionMenu(mainframe, tkvar, *Options)
    Type = "Choose the file type of the " + Title
    Label(mainframe, text= Type).grid(row = 1, column = 1)
    popupMenu.grid(row = 2, column =1)

    def ok():
        root.quit()
        root.withdraw()

    Button2 = Button(root, text="OK", command=ok)
    Button2.pack()

    root.mainloop()

    Filetype = tkvar.get()
    Filetype = Filetype.lower()

    messagebox.showinfo(Title, Message)
    Filename  = askopenfilename()

    File_Info = [ Filename, Filetype ]

    root.withdraw()
    return File_Info

def Get_Sequence():
    Target_Seq_File = Get_Files("Target Sequence", False)
    Genome_Seq_File = Get_Files("Genome Sequence", False)

    Target = SeqIO.read(Target_Seq_File[0], Target_Seq_File[1] )
    Genome = SeqIO.read(Genome_Seq_File[0], Genome_Seq_File[1] )

    Additional_Seq_File = Get_Files("Additional Sequence", True)
    if not(Additional_Seq_File == 1):
        Genome = Genome + SeqIO.read(Additional_Seq_File[0], Additional_Seq_File[1] )

    while not (Additional_Seq_File == 1 ):
         Additional_Seq_File = Get_Files("Additional Sequence", True)
         if not(Additional_Seq_File == 1):
            Genome = Genome + SeqIO.read(Additional_Seq_File[0], Additional_Seq_File[1] )

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

messagebox.showinfo("Searching", "Please Wait")
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
