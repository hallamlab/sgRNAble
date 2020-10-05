"""
Module for generating list of potentail guides
"""
import numpy as np
from Bio.Seq import Seq
from optimal_guide_finder.Azimuth_Model import model_comparison

PAM = "GG"
GUIDE_RNA_LENGTH = 20
AZIMUTH_DISTANCE = 24
TOTAL_AZIMUTH_DISTANCE = 30

def find_pams(sequence):
    """
    Find the Guide RNAs in a Sequence
    """
    locations = []
    sequence = str(sequence)
    position = 0

    while position < len(sequence):
        i = sequence[position:].find(PAM)
        if i < 0:
            break

        potential_guide_location = position + i - GUIDE_RNA_LENGTH

        # check if the guide is too close to the beginning of the gene.
        # Plus one is due to the N in the GG so the entire frame should be shifted down
        if potential_guide_location < (AZIMUTH_DISTANCE - GUIDE_RNA_LENGTH) + 1:
            pass

        # check if the guide is long enough for azimuth analysis
        elif potential_guide_location + TOTAL_AZIMUTH_DISTANCE < len(sequence):
            # the negative 1 accounts for the N in the GG
            locations.append(potential_guide_location - 1)

        position = position + i + 1

    return locations

def select_guides(target_dict, args):
    """
    TODO: Seperate the different elif statements into functions
    TODO: Don't pass in entire args object, pass in fields needed
    """
    # Obtain the Guide RNAs from the Target Sequence
    guide_list = {}

    # Default purpose, the tool looks for all avalible guides possibly found in the sequence given by the user
    if args.purpose == "d":
        for gene in target_dict:
            # Get guides on the positive strand
            guide_list[gene] = []
            locations = find_pams(target_dict[gene])
            strand_array = ["Positive"] * len(locations)

            # Get guides on the negative strand
            neg_locations = find_pams(Seq.reverse_complement(target_dict[gene]))
            sequence_length = len(target_dict[gene])
            neg_locations = [sequence_length - x for x in neg_locations]
            negative_array = ["Negative"] * len(neg_locations)

            # Combine the information into a single array for each gene
            # Magic numbers increase the sequence by 6 bps after NGG and 4 to raise total to 30
            guides = [str(target_dict[gene][loc-4:loc+GUIDE_RNA_LENGTH+6]) for loc in locations]

            # Magic numbers perform same as above, but in the opposite direction due to negative strand
            guides.extend([str(Seq.reverse_complement(
                target_dict[gene][loc-6-GUIDE_RNA_LENGTH:loc+4])) for loc in neg_locations])
            locations.extend(neg_locations)
            strand_array.extend(negative_array)

            #In case the gene is too short to have no Pam sites
            if not guides:
                guide_list.pop(gene)
                continue

            # Run the model through the Azimuth Model
            predictions = model_comparison.predict(np.array(guides))
            guides = [x for y, x in sorted(zip(predictions, guides), reverse=True)][:args.azimuth_cutoff]
            locations = [x for y, x in sorted(zip(predictions, locations), reverse=True)][:args.azimuth_cutoff]
            strand_array = [x for y, x in sorted(zip(predictions, strand_array), reverse=True)][:args.azimuth_cutoff]
            guides = [guide[4:24] for guide in guides]
            guide_list[gene].extend([guides, locations, strand_array])

    # only runs the negative strand as CRISPRi works better on negative strands
    elif args.purpose == "i":
        for gene in target_dict:
            # Get guides on the negative strand
            guide_list[gene] = []
            locations = find_pams(Seq.reverse_complement(target_dict[gene]))
            sequence_length = len(target_dict[gene])
            locations = [sequence_length - x for x in locations]
            strand_array = ["Negative"] * len(locations)

            # Combine the information into a single array for each gene
            # Magic numbers increase the sequence by 6 bps after NGG and 4 to raise total to 30
            guides = ([str(Seq.reverse_complement(
                target_dict[gene][loc-6-GUIDE_RNA_LENGTH:loc+4])) for loc in locations])

            #In case the gene is too short to have no Pam sites
            if not guides:
                guide_list.pop(gene)
                continue

            # Run the model through the Azimuth Model
            predictions = model_comparison.predict(np.array(guides))
            guides = [x for y, x in sorted(zip(predictions, guides), reverse=True)][:args.azimuth_cutoff]
            locations = [x for y, x in sorted(zip(predictions, locations), reverse=True)][:args.azimuth_cutoff]
            strand_array = [x for y, x in sorted(zip(predictions, strand_array), reverse=True)][:args.azimuth_cutoff]
            guides = [guide[4:24] for guide in guides]
            guide_list[gene].extend([guides, locations, strand_array])

    elif args.purpose == "g":
        pass

    return guide_list
