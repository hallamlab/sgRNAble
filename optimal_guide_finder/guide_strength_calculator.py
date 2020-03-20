import math
from multiprocessing import Process, Queue
import numpy as np
import pandas as pd
from cas_model import CasModel

NT_POS = {'A':0, 'T':1, 'C':2, 'G':3}

def initalize_model(guide_info, filename):
    """
    return a pandas dataframe with all the data:
        - same order as dictionary + source and target position in the beginning
    """
    model = CasModel(filename)
    num_offsite_targets = 0

    for (source, _) in model.genome_dictionary.items():
        for full_pam in model.get_all_pams():
            num_offsite_targets += len(model.genome_dictionary[source][full_pam])

    print("num offsite targets\n", num_offsite_targets)

    q = Queue()

    threads = []
    for gene in guide_info:
        for i, guide in enumerate(guide_info[gene][0]):
            # call
            process = Process(target=process_guide, args=(model, guide, i, q))
            process.start()
            threads.append(process)

    # pull results from the queue
    frames = []
    for gene in guide_info:
        for i, guide in enumerate(guide_info[gene][0]):
            frames.append(pd.DataFrame.from_records(q.get()))

    df = pd.concat(frames)

    # df.to_csv('test.csv', index=False)

    return model, df

def process_guide(model, guide, guide_index, queue):
    print(guide)

    num_guide = np.array([NT_POS[nt] for nt in list(guide)])
    partition_function = 1

    res = {}

    result = []

    for (source, _) in model.genome_dictionary.items():
        res[source] = {}

        for full_pam in model.get_all_pams():
            dg_pam = model.calc_dg_pam(full_pam)
            dg_supercoiling = model.calc_dg_supercoiling(sigma_initial=-0.05, target_seq=20 * "N")

            for (target_sequence, target_position) in model.genome_dictionary[source][full_pam]:
                dg_exchange = model.calc_dg_exchange(num_guide, target_sequence)
                dg_target = dg_pam + dg_supercoiling + dg_exchange

                result.append([source, target_position, target_sequence, dg_pam, full_pam,
                               dg_exchange, dg_supercoiling, dg_target])
                partition_function += math.exp(-dg_target / model.RT)

    print('\t' + "No." + str(guide_index + 1))
    print('\t' + guide)

    queue.put(result)
    return

def select_guides(model, guide_data, verbose=False):
    return []
