from Cas9_Calculator_py2 import *

if __name__ == "__main__":

    guideSequence = 'TACGTACACAAGAGCTCTAG'
    Cas9Calculator=clCas9Calculator(['/Users/siddarthraghuvanshi/Documents/GitHub/guideRNA_finder/Tiplasmidsequence.fasta'])
    sgRNA1 = sgRNA(guideSequence, Cas9Calculator)
    sgRNA1.run()
    sgRNA1.exportAsDill()
    sgRNA1.printTopTargets()
