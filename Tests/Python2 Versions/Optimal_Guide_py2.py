from __future__ import absolute_import
import numpy
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import Seq
from Tkinter import *
from tkFileDialog import askopenfilename
import tkMessageBox
from Cas9_Calculator_py2 import *


def Get_Files (Title, Repeat):

    root = Tk()
    Message = u"Please select the " + Title + u" file"
    More_Seqs = True

    # Add a grid
    mainframe = Frame(root)
    mainframe.grid(column=0,row=0, sticky=(N,W,E,S) )
    mainframe.columnconfigure(0, weight = 1)
    mainframe.rowconfigure(0, weight = 1)
    mainframe.pack(pady = 100, padx = 100)

    if(Repeat):
        Message = u"Please select an " + Title + u" file"
        More_Seqs = tkMessageBox.askyesno( Title,u"Would you like to add Additional Sequences to the genome \n (Such as Plasmids)")
        if not(More_Seqs):
            return 1

    tkMessageBox.showinfo(Title, u"Please select the file type")

    # Create a Tkinter variable
    tkvar = StringVar(root)

    # Dictionary with options
    Options = [ u"Fasta", u"Genebank",]
    tkvar.set(u'Fasta') # set the default option

    popupMenu = OptionMenu(mainframe, tkvar, *Options)
    Type = u"Choose the file type of the " + Title
    Label(mainframe, text= Type).grid(row = 1, column = 1)
    popupMenu.grid(row = 2, column =1)

    def ok():
        root.quit()
        root.withdraw()

    Button2 = Button(root, text=u"OK", command=ok)
    Button2.pack()

    root.mainloop()

    Filetype = tkvar.get()
    Filetype = Filetype.lower()

    tkMessageBox.showinfo(Title, Message)
    Filename  = askopenfilename()

    File_Info = [ Filename, Filetype ]

    root.withdraw()
    return File_Info

def Get_Sequence():
    Target_Seq_File = Get_Files(u"Target Sequence", False)
    Genome_Seq_File = Get_Files(u"Genome Sequence", False)

    Target = SeqIO.read(Target_Seq_File[0], Target_Seq_File[1] )
    Genome = SeqIO.read(Genome_Seq_File[0], Genome_Seq_File[1] )

    Additional_Seq_File = Get_Files(u"Additional Sequence", True)
    if not(Additional_Seq_File == 1):
        Genome + SeqIO.read(Additional_Seq_File[0], Additional_Seq_File[1] )

    while not (Additional_Seq_File == 1 ):
         Additional_Seq_File = Get_Files(u"Additional Sequence", True)
         if not(Additional_Seq_File == 1):
             Genome + SeqIO.read(Additional_Seq_File[0], Additional_Seq_File[1] )

         return Target.seq.upper(), Genome

#Length of the Guide RNA desired
Guide_RNA_length = 20

#Find the Guide RNAs in a Sequence
def PAM_Finder(Sequence, PAM, Direction):
  Guide_RNAs = []

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
            Guide_RNAs.append(Sequence[Position-23:Position-3])
        if(Direction < 0):
            Guide_RNAs.append(Sequence[Position+1:Position+21])
    Temp_Sequence = Temp_Sequence[i+2:]
    j = j+1

  return Guide_RNAs

#Combine the Coding and Template Strands into a single strand
def CombinetoStr (Template_Guides, Coding_Guides):
  Guides = []

  for i in xrange (len(Template_Guides)):
    if (i < len(Template_Guides)):
      Guides.append(unicode(Template_Guides[i]))

  for i in xrange (len(Coding_Guides)):
    if (i < len(Coding_Guides)):
      Guides.append(ReverseComplement(unicode(Coding_Guides[i])))

  return Guides

def ReverseComplement(nucleotide_sequence):
  comp = []
  for c in nucleotide_sequence:
    if c == u'A' or c == u'a':
      comp.append(u'T')
    if c == u'G' or c == u'g':
      comp.append(u'C')
    if c == u'U' or c == u'u' or c == u'T' or c == u't':
      comp.append(u'A')
    if c == u'C' or c == u'c':
      comp.append(u'G')
  rev_comp = u''.join(reversed(comp))
  return rev_comp

Target_Seq, Genome = Get_Sequence()

Genome = Genome + Genome.reverse_complement()
SeqIO.write(Genome, u"Total_Genome_Plus_RC", u"fasta")

tkMessageBox.showinfo(u"Searching", u"Please Wait")
#Obtain the Guide RNAs from the Target Sequence
T_Guides_GG = PAM_Finder(Target_Seq, u"GG",1)
T_Guides_CC = PAM_Finder(Target_Seq, u"CC", -1)

Target_Guides = CombinetoStr(T_Guides_GG, T_Guides_CC)

i = 0

for Guide in Target_Guides:

        print Guide
        Cas9Calculator=clCas9Calculator([u'Total_Genome_Plus_RC'])
        sgRNA1 = sgRNA(Guide, Cas9Calculator)
        sgRNA1.run()
        sgRNA1.printTopTargets()
        i += 1
        if (i == 2 ):
            break
