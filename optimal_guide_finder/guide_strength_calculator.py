import math
import numpy as np
from cas_model import CasModel

NT_POS = {'A':0, 'T':1, 'C':2, 'G':3}

MTYPE = 'U25, i4, U25, i4, U25, f8, f8, f8, i4'

def initalize_model(guide_info, filename):
    """
    return a 2D array saved in class variable:
        - same order as dictionary + source and target position in the beginning
    """
    model = CasModel(filename)
    num_offsite_targets = 0

    for (source, _) in model.genome_dictionary.items():
        for full_pam in model.get_all_pams():
            num_offsite_targets += len(model.genome_dictionary[source][full_pam])

    print("num offsite targets\n", num_offsite_targets)

    for gene in guide_info:
        calculate_binding(model, guide_info, gene)

    return model

def calculate_binding(model, guide_info, gene, verbose=True):
    for i, guide in enumerate(guide_info[gene][0]):
        print(guide)

        num_guide = np.array([NT_POS[nt] for nt in list(guide)])
        target_sequence_energetics = np.zeros(shape=(len(model.genome_dictionary.items()), ), dtype=MTYPE)

        partition_function = 1

        for j, (source, _) in enumerate(model.genome_dictionary.items()):
            for full_pam in model.get_all_pams():
                dg_pam = model.calc_dg_pam(full_pam)
                dg_supercoiling = model.calc_dg_supercoiling(sigma_initial=-0.05, target_seq=20 * "N")

                for (target_sequence, target_position) in model.genome_dictionary[source][full_pam]:
                    dg_exchange = model.calc_dg_exchange(num_guide, target_sequence)
                    dg_target = dg_pam + dg_supercoiling + dg_exchange

                    target_sequence_energetics[j] = np.array([(
                        source, target_position, target_sequence, dg_pam, full_pam,
                        dg_exchange, dg_supercoiling, dg_target, partition_function
                    )], dtype=MTYPE)

                    partition_function += math.exp(-dg_target / model.RT)

        if verbose:
            print('\t' + "No." + str(i + 1))
            print('\t' + guide)
            print('\t' + "Position in Target Seq:" + str(guide_info[gene][1][i]))
            print('\t' + "Strand: " + str(guide_info[gene][2][i]) + '\n')

def select_guides(model, guide_data, verbose=False):
    return []
