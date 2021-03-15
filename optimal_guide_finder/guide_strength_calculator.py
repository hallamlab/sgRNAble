"""
Module for calculating off target binding effects
"""
import math
import time
import logging
from multiprocessing.pool import ThreadPool
import numpy as np
import pandas as pd
from optimal_guide_finder.cas_model import CasModel

NT_POS = {'A':0, 'T':1, 'C':2, 'G':3}

def initalize_model(filename):
    """
    Init Cas9 binding model

    Args:
        filename (str): Path to genome file generated using get_sequence

    Returns:
        CasModel
    """
    logger = logging.getLogger(__name__)

    __start = time.time()
    logger.info("Creating Model...")
    model = CasModel(filename)
    __elasped = (time.time() - __start)
    logger.info("Time Model Building: {:.2f}".format(__elasped))

    return model

def process_guides(model, guide_info, num_threads=None):
    """
    Process guides for off target binding effects using model for Cas9 binding

    Args:
        model (CasModel): Initialized model object
        guide_info (dict): Dictionary gernerated by guide generator containing guides information
        num_threads (int, optional): max number of threads to use, None will use all available threads.
                                     Defaults to None

    Returns:
        Dataframe: Pandas dataframe containing results
    """
    logger = logging.getLogger(__name__)

    info_df = pd.DataFrame()
    logger.info("Processing Guides...")
    __start = time.time()

    pool = ThreadPool(processes=num_threads)
    results = []
    for gene in guide_info:
        for i, guide in enumerate(guide_info[gene][0]):
            guide_data = pd.Series([guide, gene, guide_info[gene][1][i], guide_info[gene][2][i]],
                                   index=["Guide Sequence", "Gene/ORF Name", "Location in Gene", "Strand"])

            # call
            res = pool.apply_async(_process_guide, (model, guide, logger))
            results.append(res)
            info_df = info_df.append(guide_data, ignore_index=True)

    pool.close()

    # pull results from the queue
    result_df = pd.DataFrame()
    for res in results:
        result_df = result_df.append(res.get(), ignore_index=True)

    results_df = pd.merge(info_df, result_df, on='Guide Sequence')

    __elasped = (time.time() - __start)
    logger.info("Time Spent Analysing Guides: {:.2f}".format(__elasped))

    return results_df

def _process_guide(model, guide, logger):
    """
    Process a single guide for off target binding effects

    Args:
        model (CasModel): initialized model
        guide (string): guide sequence
        logger (logging.logger): logger instance to use

    Returns:
        [Series]: Pandas Series containing results for guide
    """
    num_guide = np.array([NT_POS[nt] for nt in list(guide)])
    partition_function = 1

    result = []
    dgs = []

    for (source, _) in model.genome_dictionary.items():
        for full_pam in model.get_all_pams():
            dg_pam = model.calc_dg_pam(full_pam)
            dg_supercoiling = model.calc_dg_supercoiling(sigma_initial=-0.05, target_seq=20 * "N")

            for (target_sequence, _) in model.genome_dictionary[source][full_pam]:
                np_target_sequence = np.array([NT_POS[nt] for nt in list(target_sequence)])
                dg_exchange = model.calc_dg_exchange(num_guide, np_target_sequence)
                dg_target = dg_pam + dg_supercoiling + dg_exchange

                result.append([target_sequence, math.exp(-dg_target / model.RT)])
                partition_function += math.exp(-dg_target / model.RT)
                dgs.append(dg_target)

    result.insert(0, [guide, partition_function])
    guide_series = _process_off_target_guides(result)
    logger.info(guide_series)
    dgs.sort()
    info_dict = _info_logging(partition_function, dgs[:5], model.RT)
    info_dict['guide'] = guide
    logger.info(info_dict)

    return guide_series

def _process_off_target_guides(guide_result):
    """
    Post processing on guide result to calculate entropy and summarize results

    Args:
        guide_result (list): result from processing guides

    Returns:
        [Series]: Pandas series summarizing results for guide
    """
    guide_seq = guide_result[0][0]
    partition_function = guide_result[0][1]
    guide_entropy = 0
    exact_matches = 0
    for off_target in guide_result[1:]:
        if guide_seq == off_target[0]:
            exact_matches += 1
        probability = off_target[1]/partition_function
        guide_entropy -= probability*np.log2(probability)

    guide_series = pd.Series([guide_seq,
                              guide_entropy,
                              exact_matches],
                             index=["Guide Sequence",
                                    "Entropy Score",
                                    "Number of Exact Matches"])
    return guide_series

def _info_logging(partition_function, dgs, RT):
    """
    Summarizes information into dictionary for logging
    """
    info_dict = {'dg':[], 'probability':[]}

    for dg in dgs:
        info_dict['dg'].append(dg)
        info_dict['probability'].append(math.exp(-dg/RT)/partition_function)

    return info_dict
