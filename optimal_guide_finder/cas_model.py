"""
Module responsible for calculating the per guide binding
"""
import numpy as np
import scipy.io
from Bio import SeqIO
from numba import types
from numba import jit, float64, int32
import logging
logger = logging.getLogger(__name__)

class CasModel():

    RT = 0.61597

    # the PAMs with the highest dG, ignoring other PAM sequences by setting their dG to 0
    _PAM_ENERGY = {
        'GGA': -9.8, 'GGT': -10, 'GGC': -10, 'GGG': -9.9, 'CGG': -8.1, 'TGG': -7.8, 'AGG': -8.1,
        'AGC': -8.1, 'AGT': -8.1, 'AGA': -7.9, 'GCT': -7.1, 'GCG': -6.9, 'ATT': -7.2, 'ATC': -6.4,
        'TTT': -7.6, 'TTG': -6.8, 'GTA': -7.4, 'GTT': -7.9, 'GTG': -7.7, 'AAT': -7, 'AAG': -7,
        'TAT': -7.2, 'TAG': -7.2, 'GAA': -7.2, 'GAT': -7.3, 'GAC': -7.2, 'GAG': -7.3
    }

    def __init__(self, filename, quick_mode=True, model_name='data/InvitroModel.mat'):
        """
        Initialize a CasCalculator object

        Arguments:
            filename {str} -- TODO

        Keyword Arguments:
            quick_mode {bool} -- TODO (default: {True})
            model_name {str} -- path to mat file containing model data (default: {'data/InvitroModel.mat'})
        """
        self._quick_mode = quick_mode
        self._model_name = model_name
        self._nt_mismatch_in_first8_list = []

        data = scipy.io.loadmat(self._model_name)

        # TODO: Get rid of this field eventually
        self._weights = data['w1']
        self._new_weights = np.array([weight[0]*2.3 for weight in self._weights])
        self._dec_nn = data['decNN']

        self._init_genome_finder(filename)

    def print_model_info(self):
        """
        TODO: refactor this function to work with the current format of weights arr
        """
        m = 0
        s = 0
        negative_val = 0
        for i, l in enumerate(self._dec_nn):
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

    def get_all_pams(self):
        # PAM part will be 'GGT'
        for (pam_part, _) in sorted(list(self._PAM_ENERGY.items()), key=lambda x: x[1]):
            for nt in ('A', 'G', 'C', 'T'):  # nt + PAMpart will be all possible 'NGGT'
                yield nt + pam_part

    def calc_dg_pam(self, pam_full_seq):
        """
        TODO: This function is not used any where
        """
        # PAM sequence is 5' - N xxx N - 3' where the energy of xxx is listed below.
        # a normal PAM of 'NGG' with 'TC' afterwards would be listed as 'GGT'
        key = pam_full_seq[1:4]
        if key in self._PAM_ENERGY:
            return self._PAM_ENERGY[key]
        else:
            return 0.0

    def calc_dg_exchange(self, guide_seq, target_seq):
        self._nt_mismatch_in_first8_list = []

        if self._quick_mode:
            solverfunc = self._quick_calc_exchange_energy
        else:
            # solverfunc = self._calc_exchange_energy (Refer to older commits for this version of the function)
            pass

        dg_exchange = solverfunc(self._new_weights, guide_seq, target_seq)

        return dg_exchange

    def calc_dg_supercoiling(self, sigma_initial, target_seq):
        sigma_final = -0.08
        dg_supercoiling = 10.0 * \
            len(target_seq) * self.RT * (sigma_final**2 - sigma_initial**2)
        return dg_supercoiling

    def _init_genome_finder(self, filename):
        genome_dictionary = {}

        handle = open(filename, 'r')
        records = SeqIO.parse(handle, "fasta")
        record = next(records)
        handle.close()

        full_sequence = str(record.seq)

        positions_at_mers = self._identify_nucleotide_positions_of_mers(full_sequence, 10)

        genome_dictionary[filename] = {}
        target_sequence_list = []

        for full_pam in self.get_all_pams():
            target_sequence_list = \
                self._identify_target_sequences_matching_pam(full_pam, positions_at_mers, full_sequence)
            genome_dictionary[filename][full_pam] = target_sequence_list

        self.genome_dictionary = genome_dictionary

    @staticmethod
    #@jit('float64(float64[:], int32[:], int32[:])', nopython=True)
    @jit(nopython=True)
    def _quick_calc_exchange_energy(weights, cr_rna, target_seq):
        """
        calculate the delta G value a potential guide against the target sequence
        this method is made static to allow compilation with numba for speed-up

        Arguments:
            weights {np-arr} -- [description]
            cr_rna {np-arr} -- [description]
            target_seq {[type]} -- [description]

        Returns:
            [folat] -- delta G value
        """
        dg = 0
        for i in range(len(cr_rna)): 
            pos = 20 - i # TODO: magic number
            if cr_rna[i] == target_seq[i]:
                continue
            else:
                dg += weights[pos]
        return float(dg)

    def _mers(self, length):
        """
        Generates multimers for sorting through list of 10mers based on user
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
                        seq_list.append(seq + x)
            counter += 1

        last_n_mers = 4**length
        return seq_list[len(seq_list)-last_n_mers:]

    def _identify_nucleotide_positions_of_mers(self, full_sequence, length=10):
        """
        Saves list of nucleotide positions in genome that all match a unique N-mer
        sequence. Counting begins at __ending_ of MER.

        Usage:  genomePositionsAtMers[mer_sequence] is a list of nucleotide positions
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

        counter = 0

        while counter < (len(full_sequence)-length):
            word = full_sequence[counter: counter+length]
            try:
                positions_at_mers[word].append(counter+length)
            except:
                logger.error('Genome sequence contains non-nucleotide character')
            counter += 1

        return positions_at_mers

    def _identify_target_sequences_matching_pam(self, pam_seq, positions_at_mers, full_sequence,
                                                target_sequence_length=20):
        """
        Generates a list of target nucleotide sequences and corresponding nt
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
