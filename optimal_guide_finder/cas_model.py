"""
TODO:
- Reconsider class naming
- Can also use some refactoring
- Pull out some helper methods from the run function
- Revise class fields
"""
import math
from time import time
from cas_calculator import CasCalculator
import numpy as np

class sgRNA():

    _PRINT = True
    _NUM_TARGETS_RETURNED = 5

    def __init__(self, guide_info, filename):
        self.guide_info = guide_info
        self._cas_calculator = CasCalculator(filename)
        self._target_sequence_energetics = {}
        self._partition_function = 1
        self._debug = False

    def run(self):
        genome_dictionary = self._cas_calculator.genome_dictionary
        num_offsite_targets = 0
        nt_pos={'A':0,'T':1,'C':2,'G':3}

        for (source, targets) in genome_dictionary.items():
            for full_pam in self._cas_calculator.get_all_pams():
                num_offsite_targets += len(genome_dictionary[source][full_pam])

        print("num offsite targets\n", num_offsite_targets)

        for gene in self.guide_info:
            for i, guide in enumerate(self.guide_info[gene][0]):
                print(guide)
                num_guide = np.array([nt_pos[nt] for nt in list(guide)])

                begin_time = time()

                self._partition_function = 1

                for (source, targets) in genome_dictionary.items():
                    self._target_sequence_energetics[source] = {}
                    for full_pam in self._cas_calculator.get_all_pams():
                        dg_pam = self._cas_calculator.calc_dg_pam(full_pam)
                        #dG_PAM = 0
                        #dG_supercoiling= 0
                        dg_supercoiling = self._cas_calculator.calc_dG_supercoiling(
                            sigmaInitial=-0.05, targetSequence=20 * "N")  # only cares about length of sequence
                        for (target_sequence, target_position) in genome_dictionary[source][full_pam]:
                            dg_exchange = self._cas_calculator.calc_dG_exchange(num_guide, target_sequence)
                            #dG_exchange = 0
                            dg_target = dg_pam + dg_supercoiling + dg_exchange

                            self._target_sequence_energetics[source][target_position] = {
                                'sequence': target_sequence,
                                'dG_PAM': dg_pam,
                                'full_PAM': full_pam,
                                'dG_exchange': dg_exchange,
                                'dG_supercoiling': dg_supercoiling,
                                'dG_target': dg_target
                            }

                            self._partition_function += math.exp(-dg_target / self._cas_calculator.RT)

                if self._PRINT:
                    print('\t' + "No." + str(i + 1))
                    print('\t' + guide)
                    print('\t' + "Position in Target Seq:" +
                          str(self.guide_info[gene][1][i]))
                    print('\t' + "Strand: " +
                          str(self.guide_info[gene][2][i]) + '\n')
                # print("\n")
                #print('\t'.join( [	"Position in Genome", "Binding site", "dG_Target", "Partition Function"] ))

                for (source, targets) in list(self._target_sequence_energetics.items()):
                    #print("SOURCE: %s" % source)

                    # sort by smallest to largest dG_target
                    sorted_target_list = sorted(list(targets.items()), key=lambda k_v: k_v[1]['dG_target'])

                    if self._PRINT:
                        print("POSITION\tTarget Sequence\tdG_Target\t% Partition Function")

                    for (position, info) in sorted_target_list[0:self._NUM_TARGETS_RETURNED]:
                        percent_partition_function = 100 * math.exp(-info['dG_target'] / self._cas_calculator.RT) / self._partition_function
                        if self._PRINT:
                            print("%s\t%s\t%s\t%s" % (str(position),
                                                      (" "*3 +
                                                       info['sequence']),
                                                      str(round(
                                                          info['dG_target'], 2)),
                                                      str(percent_partition_function)))

                end_time = time()

                print("Elapsed Time: {:.2f}".format(end_time - begin_time))
                # print()
                if self._PRINT:
                    print("\n\n")
