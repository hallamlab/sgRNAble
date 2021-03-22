"""
Module responsible for calculating the per guide binding
"""
import os
import logging
import numpy as np
import scipy.io
from Bio import SeqIO
from numba import jit
from tqdm import tqdm

class CasModel():
    """
    Model for Cas9 binding based on biophysical properties of DNA binding
    """
    RT = 0.61597

    # the PAMs with the highest dG, ignoring other PAM sequences by setting their dG to 0
    _PAM_ENERGY = {
        'GGA': -9.8, 'GGT': -10, 'GGC': -10, 'GGG': -9.9, 'CGG': -8.1, 'TGG': -7.8, 'AGG': -8.1,
        'AGC': -8.1, 'AGT': -8.1, 'AGA': -7.9, 'GCT': -7.1, 'GCG': -6.9, 'ATT': -7.2, 'ATC': -6.4,
        'TTT': -7.6, 'TTG': -6.8, 'GTA': -7.4, 'GTT': -7.9, 'GTG': -7.7, 'AAT': -7, 'AAG': -7,
        'TAT': -7.2, 'TAG': -7.2, 'GAA': -7.2, 'GAT': -7.3, 'GAC': -7.2, 'GAG': -7.3
    }

    def __init__(self, filename, model_name='data/InvitroModel.mat'):
        """
        Initialize a CasModel object

        Arguments:
            filename {string} -- location of saved genome and target

        Keyword Arguments:
            model_name {str} -- path to mat file containing model data (default: {'data/InvitroModel.mat'})
        """
        root_path = os.path.dirname(__file__)
        model_name = os.path.join(root_path, model_name)
        self._logger = logging.getLogger(__name__)
        self._model_name = model_name
        self._nt_mismatch_in_first8_list = []

        data = scipy.io.loadmat(self._model_name)

        self._weights = data['w1']
        self._new_weights = np.array([weight[0]*2.3 for weight in self._weights])
        self._dec_nn = data['decNN']

        self._init_genome_finder(filename)

    def print_model_info(self):
        """
        Prints current model information
        """
        count = 0
        sums = 0
        negative_val = 0
        for i, mismatch in enumerate(self._dec_nn):
            for j, energy in enumerate(mismatch):
                if float(energy) < 0:
                    negative_val += 1
                if i != j:
                    sums += float(energy)
                    count += 1
        mean_nn = float(sums)/float(count)

        sum_weights = 0
        for weight in self._weights:
            sum_weights += weight
        meanw = sum_weights/len(self._weights)

        self._logger.info("Average mismatch energy: %d", mean_nn)
        self._logger.info("Average weight: %d", meanw)
        self._logger.info("Number of negative energies: %d", negative_val)

    def get_all_pams(self):
        """
        Generates a list of all possible PAMs

        Yields:
            [string]: Pam seqs
        """
        # PAM part will be 'GGT'
        for (pam_part, _) in sorted(list(self._PAM_ENERGY.items()), key=lambda x: x[1]):
            for nucleotide in ('A', 'G', 'C', 'T'):  # nt + PAMpart will be all possible 'NGGT'
                yield nucleotide + pam_part

    def calc_dg_pam(self, pam_full_seq):
        """
        Calculates the delta-G for full pam sequence
        """
        # PAM sequence is 5' - N xxx N - 3' where the energy of xxx is listed below.
        # a normal PAM of 'NGG' with 'TC' afterwards would be listed as 'GGT'
        key = pam_full_seq[1:4]
        if key in self._PAM_ENERGY:
            return self._PAM_ENERGY[key]

        return 0.0

    def calc_dg_exchange(self, guide_seq, target_seq):
        """
        Calculate the Delta-G Exchange for guide and target sequences
        """
        self._nt_mismatch_in_first8_list = []

        solverfunc = self._quick_calc_exchange_energy
        dg_exchange = solverfunc(self._new_weights, guide_seq, target_seq)

        return dg_exchange

    def calc_dg_supercoiling(self, sigma_initial, target_seq):
        """
        Calculate the Delta-G supercoiling for guide and target sequences
        """
        sigma_final = -0.08
        dg_supercoiling = 10.0 * \
            len(target_seq) * self.RT * (sigma_final**2 - sigma_initial**2)
        return dg_supercoiling

    def _init_genome_finder(self, filename):
        """
        Initializes genome dictionary by finding all potrntial locations in genome for off-target binding

        Args:
            filename (string): Path to genome and target files
        """
        genome_dictionary = {}

        handle = open(filename, 'r')
        records = SeqIO.parse(handle, "fasta")
        record = next(records)
        handle.close()

        full_sequence = str(record.seq)

        positions_at_mers = self._identify_nucleotide_positions_of_mers(full_sequence, 10)

        self._logger.info("Identifying target sites...")
        genome_dictionary[filename] = {}
        target_sequence_list = []
        mers_omitted_tally = 0

        for full_pam in tqdm(self.get_all_pams()):
            target_list_output = \
                self._identify_target_sequences_matching_pam(full_pam, positions_at_mers, full_sequence)
            target_sequence_list = target_list_output[0]
            mers_omitted_tally += target_list_output[1]
            genome_dictionary[filename][full_pam] = target_sequence_list

        if mers_omitted_tally > 0:
            self._logger.error('%d k-mers with non-nucleotide characters omitted from target list', mers_omitted_tally)

        self.genome_dictionary = genome_dictionary

    @staticmethod
    @jit(nopython=True)
    def _quick_calc_exchange_energy(weights, cr_rna, target_seq):
        """
        Calculate the delta G value a potential guide against the target sequence
        this method is made static to allow compilation with numba for speed-up
        """
        dg = 0
        for i in range(len(cr_rna)):
            pos = 20 - i
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

        non_nucleotide_count = 0 # Added
        while counter < (len(full_sequence)-length):
            word = full_sequence[counter: counter+length]
            try:
                positions_at_mers[word].append(counter+length)
            except KeyError:
                non_nucleotide_count += 1
            counter += 1

        if non_nucleotide_count > 0:
            self._logger.error('%d k-mers with non-nucleotide characters identified', non_nucleotide_count)

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
        mers_omitted_count = 0
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
                    if set(target_sequence).issubset('ATCG'):
                        target_sequence_list.append((target_sequence, nt))
                        # print(target_sequence, ' : ', nt)
                    else:
                        mers_omitted_count += 1
        # print(mers_omitted_count)
        target_list_output = [target_sequence_list, mers_omitted_count]
        return target_list_output
