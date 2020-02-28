"""
TODO: This class could really use some refactoring
"""

import numpy as np
import scipy.io
from Bio import SeqIO
from numba import jit

class CasCalculator():

    RT = 0.61597

    # the PAMs with the highest dG, ignoring other PAM sequences by setting their dG to 0
    _PAM_ENERGY = {
        'GGA': -9.8, 'GGT': -10, 'GGC': -10, 'GGG': -9.9, 'CGG': -8.1, 'TGG': -7.8, 'AGG': -8.1,
        'AGC': -8.1, 'AGT': -8.1, 'AGA': -7.9, 'GCT': -7.1, 'GCG': -6.9, 'ATT': -7.2, 'ATC': -6.4,
        'TTT': -7.6, 'TTG': -6.8, 'GTA': -7.4, 'GTT': -7.9, 'GTG': -7.7, 'AAT': -7, 'AAG': -7,
        'TAT': -7.2, 'TAG': -7.2, 'GAA': -7.2, 'GAT': -7.3, 'GAC': -7.2, 'GAG': -7.3
    }

    def __init__(self, filename, quick_mode=True, model_name='data/InvitroModel.mat'):
        self._quick_mode = quick_mode
        self._model_name = model_name

        data = scipy.io.loadmat(self._model_name)
        self._weights = data['w1']
        self._new_weights = np.array([weight[0]*2.3 for weight in self._weights])
        self._decNN = data['decNN']

        self._init_genome_finder(filename)

    def get_all_pams(self):
        # PAMpart will be 'GGT'
        for (pam_part, _) in sorted(list(self._PAM_ENERGY.items()), key=lambda x: x[1]):
            for nt in ('A', 'G', 'C', 'T'):  # nt + PAMpart will be all possible 'NGGT'
                yield nt + pam_part

    def print_model_info(self):
        m = 0
        s = 0
        negative_val = 0
        for i, l in enumerate(self._decNN):
            for j, e in enumerate(l):
                if float(e) < 0:
                    negative_val += 1
                if i != j:
                    s += float(e)
                    m += 1
        meanNN = float(s)/float(m)

        sw = 0
        for w in self._weights:
            sw += w

        meanw = sw/len(self._weights)
        print('average mismatchc energy: ', meanNN)
        print('average weight:', meanw)
        print('number of negative energies: ', negative_val)

    def calc_dg_pam(self, pam_full_seq):
        # PAM sequence is 5' - N xxx N - 3' where the energy of xxx is listed below. A normal PAM of 'NGG' with 'TC' afterwards would be listed as 'GGT'
        key = pam_full_seq[1:4]
        if key in self._PAM_ENERGY:
            return self._PAM_ENERGY[key]
        else:
            return 0.0

        # acceptedPAMList=PAM_dic_energy.keys()
        # self.dG_PAM_List=[]
        # self.WarningPAM_List=[]
        # PAMsize=len(self.PAM)
        # for target in self.sequence_list:
            # tPAM=target[-(PAMsize):-1]+target[-1]
            # if tPAM in acceptedPAMList:
            # dGPAM=PAM_dic_energy[tPAM]
            # warning=''
            # else:
            # dGPAM=0
            # warning='N.B'
            # self.dG_PAM_List.append(dGPAM)
            # self.WarningPAM_List.append(warning)

    def calc_dG_exchange(self, guideSequence, targetSequence):
        self.nt_mismatch_in_first8_list = []
        if self._quick_mode:
            solverfunc = self._quick_calc_exchange_energy
        else:
            solverfunc = self._calc_exchange_energy

        dG_exchange = solverfunc(self._new_weights, guideSequence, targetSequence)

        return dG_exchange

    def calc_dG_supercoiling(self, sigmaInitial, targetSequence):

        sigmaFinal = -0.08
        dG_supercoiling = 10.0 * \
            len(targetSequence) * self.RT * (sigmaFinal**2 - sigmaInitial**2)
        return dG_supercoiling

    def _init_genome_finder(self, filename):
        genome_dictionary = {}

        handle = open(filename, 'r')
        records = SeqIO.parse(handle, "fasta")
        record = next(records)
        handle.close()

        full_sequence = str(record.seq)
        print("full seq length", len(record.seq))
        positions_at_mers = self._identify_nucleotide_positions_of_mers(full_sequence, 10)
        print("computed positions")

        genome_dictionary[filename] = {}
        target_sequence_list = []

        for full_pam in self.get_all_pams():
            print("Full PAM", full_pam)
            target_sequence_list = self._identify_target_sequences_matching_pam(full_pam, positions_at_mers, full_sequence)
            genome_dictionary[filename][full_pam] = target_sequence_list

        self.genome_dictionary = genome_dictionary

    def _calc_exchange_energy(self, crRNA, targetSeq):
        nt_pos = {'A': 0, 'T': 1, 'C': 2, 'G': 3,
                  'a': 0, 't': 1, 'c': 2, 'g': 3}
        dG = 0
        RNA = ''
        DNA = ''
        for i in range(0, len(crRNA)):
            if i > 0:
                RNA = crRNA[(i-1):(i+1)]
                DNA = targetSeq[(i-1):(i+1)]
                RNA_index = nt_pos[RNA[0]]+4*nt_pos[RNA[1]]
                DNA_index = nt_pos[DNA[0]]+4*nt_pos[DNA[1]]

                dG1 = float(self._decNN[RNA_index][DNA_index])
                if abs(dG1-0.000015) < 1e-6:
                    dG1 = 10000
                    # during model identification, I set the value of every unknown dG to 0.000015 (if I did not find a value for it)
                    dG1 = 2.3

                pos = 20-i
                w1 = float(self._weights[pos])
                # print 'b1',RNA[0],RNA[1],DNA[0],DNA[1],RNA_index, DNA_index, pos,dG1, w1
            else:
                w1 = 0
                dG1 = 0
            if i < (len(crRNA)-1):
                RNA2 = crRNA[i:(i+2)]
                DNA2 = targetSeq[i:(i+2)]
                RNA_index = nt_pos[RNA2[0]]+4*nt_pos[RNA2[1]]
                DNA_index = nt_pos[DNA2[0]]+4*nt_pos[DNA2[1]]
                dG2 = float(self._decNN[RNA_index][DNA_index])
                if abs(dG2-0.000015) < 1e-6:
                    dG2 = 10000
                    # during model identification, I set the value of every unknown dG to 0.000015
                    # (if I did not find a value for it)
                    dG2 = 2.3

                pos = 20-i-1
                w2 = float(self._weights[pos])
                # print 'b2',RNA2[0],RNA2[1],DNA2[0],DNA2[1],RNA_index, DNA_index, pos,dG2, w2
            else:
                w2 = 0
                dG2 = 0
            dG += w1*dG1+w2*dG2
        return float(dG)

    @staticmethod
    @jit(nopython=True)
    def _quick_calc_exchange_energy(weights, crRNA, target_seq):
        dG = 0
        for i in range(len(crRNA)):
            pos = 20 - i
            if crRNA[i] == target_seq[i]:
                continue
            else:
                dG += weights[pos]
        return float(dG)

    def _mers(self, length):
        """Generates multimers for sorting through list of 10mers based on user
        specification. Multimers generated act as the keys for generating a
        hashtable to eliminate undesired sequence patterns from those 10mers not
        found in the genome.

        Usage: mers(N) = 4^(N) unique Nmers
        """
        seq_list = ['']
        counter = 0
        while counter < length:
            for seq in seq_list:
                if len(seq) == counter:
                    for x in ['A', 'T', 'C', 'G']:
                        seq_list.append(seq+x)
            counter += 1
        last_N_Mers = 4**length
        return seq_list[len(seq_list)-last_N_Mers:]

    def _identify_nucleotide_positions_of_mers(self, full_sequence, length=10):
        """Saves list of nucleotide positions in genome that all match a unique N-mer
        sequence. Counting begins at __ending_ of MER.

        Usage:   genomePositionsAtMers[mer_sequence] is a list of nucleotide positions
        within the inputted fullSequence that match mer_sequence
                        If mer_sequence ends in a PAM site, then this can be used to match
                        the first N-3 nt of a guide strand plus a PAM site sequence.
        """

        # Create a list of all possible N-mers
        all_possible_mers = self._mers(length)

        # Search through the genome and add nucleotide positions for match to an N-mer
        positions_at_mers = {}
        for mer in all_possible_mers:
            positions_at_mers[mer] = []

        print("Number of Mers: ", len(list(positions_at_mers.keys())))
        counter = 0
        while counter < (len(full_sequence)-length):
            word = full_sequence[counter: counter+length]
            positions_at_mers[word].append(counter+length)
            counter += 1
        return positions_at_mers

    def _identify_target_sequences_matching_pam(self, pam_seq, positions_at_mers, full_sequence,
                                                target_sequence_length=20):
        """ Generates a list of target nucleotide sequences and corresponding nt
        positions for an inputted sequence that matched the PAM_seq.
                Uses the positionsAtMers dictionary to accelerate the identification.
                Good for large genomes.

        Usage:  listOfTargets = identifyTargetSequencesMatchingPAM('CGG',
                                                                positionsAtMers, genome_sequence)
        """
        target_sequence_list = []
        all_mers = list(positions_at_mers.keys())
        mer_length = len(all_mers[0])
        list_of_mers_with_pam = [
            mer + pam_seq for mer in self._mers(mer_length - len(pam_seq))]
        for mer_with_pam in list_of_mers_with_pam:
            nt_list = positions_at_mers[mer_with_pam]
            for nt in nt_list:
                begin = nt-target_sequence_length - len(pam_seq)
                end = nt - len(pam_seq)
                # Does not account for circular DNAs
                if begin > 0 and end < len(full_sequence):
                    target_sequence = full_sequence[begin: end]
                    target_sequence_list.append((target_sequence, nt))
        return target_sequence_list
