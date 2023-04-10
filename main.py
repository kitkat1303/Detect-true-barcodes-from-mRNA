from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
from Bio import SeqIO
from collections import Counter
import numpy as np
import sys
import numpy
import collections
numpy.set_printoptions(threshold=sys.maxsize)


# This program provides utility functions to process a fastq file and determine
# the true barcodes based on either the number of cells sequenced or
# a Gaussian density probability calculation of all sequences within a fastq file
# the output can either be the true barcodes found (whitelist) or the true barcodes
# along with the barcodes that are one hamming distance away from whitelisted barcodes


# fastq_to_sequences takes in a file name and returns a list of sequences
# if the file does not open, an error message is printed
def fastq_to_sequences(file_name):
    # store each line of the fastq file into sequences
    try:
        np.sequences = SeqIO.parse(open(file_name, "r"), "fastq")
    except IOError:
        print("The file entered is incorrect, please ensure it is in fastq format.")
        return
    return np.sequences


# get_barcodes takes in the sequences read in from the file along with a
# known anchor sequence and returns a list of barcodes that contain the anchor
# sequence
def get_barcodes(sequences, anchor_sequence, q_score):
    total_sequences = 0
    sequences_with_anchor = 0
    quality_sequences = 0
    np.all_barcodes_list = []
    if 0 < q_score < 31:
        print("Ensuring all bases have a quality score above: " + str(q_score))
    for RNA_seq in sequences:
        total_sequences += 1
        # if the sequence contains the anchor, then store the 30 characters preceding it
        try:
            index_start_anchor = str(RNA_seq).index(anchor_sequence)
            barcode_sequence = str(RNA_seq)[index_start_anchor - 30: index_start_anchor]
            if 0 < q_score < 31:
                is_quality_sequence = True
                for qualityNumber in RNA_seq.letter_annotations["phred_quality"]:
                    if qualityNumber < q_score:
                        is_quality_sequence = False
                        break
                # if it is a quality sequence then store barcode
                if is_quality_sequence is True:
                    quality_sequences += 1
                    sequences_with_anchor += 1
                    np.all_barcodes_list.append(barcode_sequence)
            # if no quality check is being performed then just add barcode
            else:
                np.all_barcodes_list.append(barcode_sequence)
                sequences_with_anchor += 1
        except ValueError:
            continue

    print("Total RNA sequences: " + str(total_sequences))
    print("Total RNA sequences with anchor: " + str(sequences_with_anchor))
    return np.all_barcodes_list


# get_barcode_count takes in the list of barcodes and
# returns a dict with the barcodes and their frequencies
def get_barcode_count(all_barcodes_list):
    np.frequency = collections.Counter(all_barcodes_list)
    return dict(np.frequency)


# get_whitelist takes in the frequency dict of all barcodes and returns the true barcodes
# based on statistical analysis using guassian density
def get_whitelist(barcode_counts, cell_num=0):
    # low abundance barcodes are filtered out (< 0.001 * the most abundant)
    threshold = 0.001 * Counter(barcode_counts).most_common(1)[0][1]
    # sort counts from highest to lowest abundance
    counts = sorted(barcode_counts.values(), reverse=True)
    # remove counts below the threshold
    counts_thresh = [x for x in counts if x > threshold]
    # convert the counts that meet the threshold criteria to log
    log_counts = np.log10(counts_thresh)

    # guassian density with hardcoded bw
    # determine kernel density estimate using Gaussian kernels
    # input is the log_counts of barcode_counts and bw_method is set to .1 as smoothing parameter
    density = gaussian_kde(log_counts, bw_method=0.1)

    # how many x values for density plot
    xx_values = 10000
    # list that stores evenly spaces values between min and max of log_counts of length 10000
    xx = np.linspace(log_counts.min(), log_counts.max(), xx_values)

    local_min = None

    # if number of cells is less than unique barcodes and greater than 0
    if 0 < cell_num < len(counts):
        # set threshold by taking top cell_num of counts
        try:
            threshold = counts[cell_num]
        except IndexError:
            print("Number of cells can not be used to determine threshold.")

    # if number of cells is not given or not possible to use to set threshold
    # use guassian density to determine rough distribution to remove outliers
    # this method is very conservative and does not produce the most accurate results
    else:
        local_mins = argrelextrema(density(xx), np.less)[0]
        local_mins_counts = []

        for pos_local_min in local_mins[::-1]:
            passing_threshold = sum([y > np.power(10, xx[pos_local_min]) for x, y in barcode_counts.items()])
            local_mins_counts.append(passing_threshold)

            if not local_min:
                if (pos_local_min >= 0.2 * xx_values and (log_counts.max() - xx[pos_local_min] > 2 or
                                                          xx[pos_local_min] < log_counts.max() / 2)):
                    local_min = pos_local_min

        # set the threshold if the number of cells could not be used
        if local_min is not None:
            threshold = np.power(10, xx[local_min])

    print("Threshold for frequency of barcodes set to: " + str(threshold))
    if cell_num or local_min is not None:
        np.final_barcodes = set([
            x for x, y in barcode_counts.items() if y > threshold])
    else:
        final_barcodes = None

    return list(np.final_barcodes)


# get_hamming_dist takes in a single true barcode from the whitelist
# along with a single barcode and returns whether they are hamming-1 away from
# one another
def get_hamming_dist(white_list_code, barcode):
    zipped_rna = zip(white_list_code, barcode)
    h_dist = 0
    for pair in zipped_rna:
        # if more than one base off, then hamming dist > 1
        if h_dist > 1:
            return h_dist
        if pair[0] != pair[1]:
            # change one base
            h_dist += 1
    return h_dist


# find_one_hamming_dist_ codes takes in a list of true barcodes (whitelist_codes)
# along with all barcodes and returns a dict that stores barcodes along with the
# whitelist codes they diverged from (hamming-1 base off)
def find_one_hamming_dist_codes(whitelist_codes, barcodes):
    np.one_hamming_dist_barcodes = {}
    for whitelist_code in whitelist_codes:
        for barcode in barcodes:
            hamming_dist = get_hamming_dist(whitelist_code, barcode)
            if hamming_dist == 1:
                np.one_hamming_dist_barcodes[barcode] = whitelist_code
    return np.one_hamming_dist_barcodes


# get_true_barcodes takes in file name in fastq format
# along with a known anchor sequence and returns a list of true
# barcodes along with a dict of barcodes that diverged from a whitelist
# barcodes by one base
def get_true_barcodes(list_sequences, anchor_sequence, expected_number_of_cells=0, get_hamming_one_off=False,
                      q_score=-1):
    all_barcodes = get_barcodes(list_sequences, anchor_sequence, q_score)
    all_barcode_counts = get_barcode_count(all_barcodes)
    true_barcodes = get_whitelist(all_barcode_counts, expected_number_of_cells)
    if get_hamming_one_off:
        one_base_off_barcodes = find_one_hamming_dist_codes(true_barcodes, all_barcodes)
        return true_barcodes, one_base_off_barcodes
    else:
        return true_barcodes


# set_analysis_cond takes in a filename, anchor sequence, and then has optional parameters to be entered
# which are number of cells that have been sequenced (set to 0 if not known), q_score to determine quality
# of bases which can be set between 0 and 30 (set to -1 if quality should not be checked), and a bool
# called get_hamming_one which if set to false, will not determine barcodes that are one hamming distance
# from true barcodes
def set_analysis_cond(filename, anchor_sequence, num_cells_sequenced=0, q_score=-1, get_hamming_one=False):
    if num_cells_sequenced < 0:
        print("Number of cells sequenced can not be negative. Will continue program with density "
              "calculation.")
        num_cells_sequenced = 0
    if q_score < 0:
        print("Quality score can not be negative. Proceeding without quality assessment.")
        q_score = -1
    if q_score > 30:
        print("Quality score can not be greater than 30. Proceeding without quality assessment.")
        q_score = -1
    return filename, anchor_sequence, num_cells_sequenced, q_score, get_hamming_one


# main
# call functions for program
# set parameters for analysis
file, anchor_seq, expect_number_cells, quality_score, get_hamming_one_off_codes = \
    set_analysis_cond("10_5_RNASeq.fastq", "GTACTGCGGCCGCTACCTA", 500, 10, False)
# store file contents
seq = fastq_to_sequences(file)
# get true barcodes along with hamming-1 barcodes from true barcodes
# the processing of this takes a bit longer due to the hamming function
# for future reference multi_processing would be implemented to increase speed
if get_hamming_one_off_codes:
    whitelist, hamming_one = get_true_barcodes(seq, anchor_seq, expect_number_cells,
                                               get_hamming_one_off_codes, quality_score)
    if expect_number_cells == 0:
        print("Number of cells sequenced: N/A")
    else:
        print("Number of cells sequenced: " + str(expect_number_cells))
    print("Number of true bar codes found: " + str(len(whitelist)))
    print("Number of hamming-1 barcodes found: " + str(len(hamming_one)))
# only get the true barcodes
else:
    whitelist = get_true_barcodes(seq, anchor_seq, expect_number_cells, get_hamming_one_off_codes, quality_score)
    if expect_number_cells == 0:
        print("Number of cells sequenced: N/A")
    else:
        print("Number of cells sequenced: " + str(expect_number_cells))
    print("Number of true bar codes found: " + str(len(whitelist)))
