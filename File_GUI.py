from tkinter import *
import numpy as np
from tkinter.filedialog import askopenfilename
import os

class GUI(object):

    def __init__(self, Master):
        self.row = 0
        self.Files = [[]]
        self.Files.append([])
        self.Filetypes = []
        self.num_of_Additional = 0
        self.Win = Master

    #Button Functions
    def End(self):
        self.Win.quit()

    def File(self,FileName, Type, Count):
        File = askopenfilename()
        FileName["text"] = os.path.basename(File)
        self.Files[0].append(File)
        if(Type == "R"):
            self.Files[1].append(Count)
        else:
            self.Files[1].append(Type)

    def Add(self):
        self.Done.grid(row = self.row + 3, column = 0, pady = 50)
        self.A_Seq.grid(row = self.row + 2, column = 0, pady = 0,)
        self.Organizer("Additional Sequnce", "R")
        self.row += 3

    def Get_Rid(self,Frame, Prompt, index):
        Frame.destroy()
        Prompt.destroy()
        self.Filetypes[index + 2] = 0

    #Function to create headings and buttons
    def Organizer(self, Heading, Type):

        #Create a Prompt and a Frame
        Prompt = Label(self.Win, text = Heading)
        F = Frame(self.Win)

        #Display the prompt
        Prompt.grid(row = self.row, column = 0,pady = 5, padx = 100)
        self.row += 1

        #IF this is an add button, follow this code to display it
        if (Type == "ADD"):
            More = Button(F,text = "Add More", command = self.Add)
            More.grid(row = self.row, column = 0)
            self.row += 1

        #IF this isn't an add button, run this code to display the buttons
        else:
            Count = self.num_of_Additional
            Formats = [ "Fasta", "Genbank"]
            Format_Var = StringVar(self.Win)
            Format_Var.set("Fasta")
            self.Filetypes.append(Format_Var)
            FileType = OptionMenu(F, Format_Var, *Formats)
            FileName = Button(F,text = "File Name", command = lambda: self.File(FileName,Type, Count))
            FileName.grid(row = 0, column = 2)
            FileType.grid(row = 0, column = 4)

            if(Type == "R"):
                Remove = Button(F,text = "Remove", command = lambda: self.Get_Rid(F,Prompt, Count))
                Remove.grid(row = 0, column = 6)
                F.grid(row = self.row, column = 0)
                self.row += 1
                self.num_of_Additional += 1


        return F

    def run(self):

        self.Win.title("Guide RNA Finder")

        self.T_Seq = self.Organizer("Please select the Target Sequence File", "TARGET")
        self.T_Seq.grid(row = self.row, column = 0)
        self.row += 1

        self.G_Seq = self.Organizer("Please select the Genome Sequence File", "GENOME")
        self.G_Seq.grid(row = self.row, column = 0)
        self.row += 1

        self.A_Seq= self.Organizer("Do you wish to add any additional sequences such as plasmids", "ADD")
        self.A_Seq.grid(row = self.row, column = 0)
        self.row += 1

        #Done Button
        self.Done = Button(self.Win, text = "Done", command = self.End)
        self.Done.grid(row = self.row, column = 0, pady = 50)
        self.row += 1


if __name__ == "__main__":

    Root = Tk()
    Program = GUI(Root)
    Program.run()
    Root.mainloop()
    print(Program.Files[1][0])
    for Format in Program.Filetypes:
        print(Format)
