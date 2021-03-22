"""
Entry point to the program

- Accept and parse user commands
- Generate Potential Guides
- Run through biophysical model and report results
"""
import os
import logging
import argparse
import time
import psutil
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from optimal_guide_finder import guide_generator, guide_strength_calculator, memory_limit

def init_parser():
    """
    Initializes a parser to get user arguments using argparse

    Returns:
        ArgumentParser -- parser ready to accept arguments
    """
    # Parser to get the files listed in the arguments
    parser = argparse.ArgumentParser(description="""This program helps you to find all possible guide RNAs that will
                                       target the gene. Then using the model created by Salis Lab,
                                       you can see the off target effects for the each possible guide.""",
                                     formatter_class=argparse.RawTextHelpFormatter)

    # Parsers to add arguements.
    parser.add_argument("-t", "--target_sequence", required=True,
                        help="The Gene Sequence of Interest (Fasta format)")
    parser.add_argument("-g", "--genome_sequence", required=True, nargs='+',
                        help="""The Genome of the organism, if targeting a plasmid, make sure to \n
                              include it as well (Fasta format)""")

    parser.add_argument("-a", "--azimuth_cutoff", type=int, required=False,
                        default=10, help="""How many guides should pass from azimuth screening,
                              the guides are passed based on descending azimuth prediction score""")
    parser.add_argument("-c", "--copy_number", required=False, default=1, nargs='+',
                        help="""Number of copies of target gene present""")
    parser.add_argument("-p", "--purpose", required=False, default="d",
                        help="""d: default, will look at all possible guides,
                            i: CRISPR interference on gene, will only look at negative strand""")

    parser.add_argument("-o", "--output_dir", required=False, default="output",
                        help="Path to output directory")

    parser.add_argument("-th", "--num_threads", required=False, default=None, type=int,
                        help="""Number of threads to use when running the program""")
    parser.add_argument("-m", "--max_memory", required=False, default=None, type=float,
                        help="""Maximum memory used by the tool in GiB,
                                defaults to using whatever memory is available in the system""")

    parser.add_argument("-v", "--verbose", action='store_true', help="""Enable verbose console logging""")

    return parser

def get_sequence(target_sequence, genome_sequences, copy_numbers):
    """
    Loads target and genome sequences from files, if multiple files are passed for genome sequences,
    they'll be combined into one file

    Args:
        target_sequence (string): path to fasta file containing target sequence
        genome_sequences (list(string)): list of paths to fasta files containing genome sequences
        copy_numbers (list(int)): number of copies of target sequence present for each target sequence

    Returns:
        [(dict, string)]: target sequence dictionary, genome sequence as string
    """
    # Reads the file using biopython and creates an object called target
    target_dict = SeqIO.to_dict(SeqIO.parse(target_sequence, "fasta"))

    for name in target_dict:
        target_dict[name] = target_dict[name].seq.upper()

    # Reads the Genome files using biopython and combines them into one genome object
    genome = SeqRecord(Seq(""))
    for i, genome_sequence in enumerate(genome_sequences):
        genome_parts = SeqIO.parse(genome_sequence, "fasta")
        for part in genome_parts:
            genome.seq = genome.seq + part.seq * int(copy_numbers[i])

    return target_dict, genome.seq.upper()

def initialize_logger(output_file, logging_level=logging.INFO):
    """
    Initialize logging to be used throughout the program

    Args:
        output_file (str): path to logging output
        logging_level (logging.LEVEL, optional): verbosity of console logger. Defaults to logging.INFO.
    """
    logging.basicConfig(level=logging.INFO,
                        filename=output_file + '/run.log',
                        filemode='w')
    console_logger = logging.StreamHandler()
    console_logger.setLevel(logging_level)
    logging.getLogger().addHandler(console_logger)

def main():
    """
    Main workflow
    """
    start_time = time.time()

    # init parser
    parser = init_parser()
    args = parser.parse_args()

    # create output directory if not already present
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    assert(os.path.isdir(args.output_dir))

    # set logger verbosity level
    if args.verbose:
        initialize_logger(args.output_dir)
    else:
        initialize_logger(args.output_dir, logging.ERROR)

    logger = logging.getLogger(__name__)

    # set memory limit if passed in
    if args.max_memory is not None:
        memory_limit.set_limit(args.max_memory)

    # create the path for genome file
    logger.info("Creating genome file...")
    genome_location = args.output_dir + '/Run_Genome'

    # Get the sequences in a Seq format from user fasta or genebank files
    logger.info("Generating target dictionary..")
    copies = [args.copy_number]*len(args.genome_sequence) if isinstance(args.copy_number, int) else args.copy_number
    target_dict, genome = get_sequence(args.target_sequence, args.genome_sequence, copies)

    ref_record = SeqRecord(genome, id="refgenome", name="reference", description="a reference background")
    ref_record = ref_record + ref_record.reverse_complement()
    SeqIO.write(ref_record, genome_location, "fasta")

    # Select the guides based on the purpose and the azimuth model
    logger.info("Selecting initial guides..")
    guide_list = guide_generator.select_guides(target_dict, args.purpose, args.azimuth_cutoff)

    # Build and run the model
    logger.info("Initializing binding strength model..")
    model = guide_strength_calculator.initalize_model(genome_location)

    logger.info("Processing guides off-target binding...")
    results_df = guide_strength_calculator.process_guides(model, guide_list, num_threads=args.num_threads)

    #generate and append Rank array
    logger.info("Generating result file..")
    results_df.sort_values(by=['Gene/ORF Name', 'Entropy Score'], inplace=True)
    results_df.drop_duplicates(inplace=True)
    rank_array = []
    for gene in results_df['Gene/ORF Name'].unique():
        num_guides = len(results_df[results_df['Gene/ORF Name'] == gene]['Guide Sequence'])
        rank_array.extend(list(np.arange(1, num_guides+1)))
    results_df['Rank in Target Gene'] = rank_array

    results_df.to_csv(args.output_dir + '/output.csv', index=False)

    process = psutil.Process(os.getpid())

    total_time = time.time() - start_time
    minutes = int(total_time//60)
    seconds = int(total_time%60)
    max_memory = process.memory_info().rss // (1024 * 1024)
    logger.info("Finished sgRNAble run, total time elapsed (min/sec): %d:%d, max memory usage (MB): %d", minutes, seconds, max_memory)

if __name__ == "__main__":
    main()
