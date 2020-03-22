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
    df = pd.DataFrame(index = ["Guide Sequence", "Entropy Score"])
    for _ in threads:
        df.append(q.get())

    #Would it be possible to pass information about the guide such as it's location in the target_seq
    #and the strand it targets to this function. That way we can apped it to each guide in this dataframe
    #df = pd.DataFrame(frames)

    df.to_csv('test.csv', index=False)

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

                result.append([math.exp(-dg_target / model.RT)])
                partition_function += math.exp(-dg_target / model.RT)
    
    results.insert(0,[guide,partition_function])
    guide_series = process_off_target_guides(results)    
    print('\t' + "No." + str(guide_index + 1))
    print('\t' + guide)

    queue.put(guide_series)
    return

def process_off_target_guides(guide_data, verbose=False):
    guide_seq = guide_info[0][0]
    partition_function = guide_info[0][1]
    guide_entropy = 0
    for off_target in guide_info[1:]:
        probability = off_target[0]/partition_function
        guide_entropy -= probability*np.log2(probability)
    guide_series = pd.Series([guide_seq,
                              guide_entropy],
                             index = ["Guide Sequence",
                                      "Entropy Score"])
    return guide_series
