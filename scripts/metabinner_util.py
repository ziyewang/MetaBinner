import os
import mimetypes
import gzip
import sys

from collections import defaultdict, Counter

import biolib.seq_io as seq_io
from Bio import SeqIO
import pandas as pd
from biolib.common import (make_sure_path_exists,
                           check_dir_exists,
                           check_file_exists,
                           query_yes_no,
                           remove_extension)


def get_length(fastx_file):
    file_type = mimetypes.guess_type(fastx_file)[1]
    if file_type == 'gzip':
        f = gzip.open(fastx_file, "rt")
    elif not file_type:
        f = open(fastx_file, "rt")
    else:
        raise RuntimeError("Unknown type of file: '{}".format(fastx_file))
    length = {}
    if os.path.getsize(fastx_file) == 0:
        return length
    file_format = None
    line = f.readline()
    if line.startswith('@'):
        file_format = "fastq"
    elif line.startswith(">"):
        file_format = "fasta"
    f.seek(0)
    if not file_format:
        raise RuntimeError("Invalid sequence file: '{}".format(fastx_file))
    for seq_record in SeqIO.parse(f, file_format):
        length[seq_record.id] = len(seq_record.seq)

    f.close()

    return length


####
####通过coverage_file得到namelist
def get_namelist(cov_file):
    namelist = pd.read_csv(cov_file, sep='\t', usecols=range(1)).values[:, 0]

    return namelist


# from checkm: from checkm.util.seqUtils import calculateN50
def calculateN50(seqLens):
    thresholdN50 = sum(seqLens) / 2.0

    seqLens.sort(reverse=True)

    testSum = 0
    N50 = 0
    for seqLen in seqLens:
        testSum += seqLen
        if testSum >= thresholdN50:
            N50 = seqLen
            break

    return N50


def cal_comp_cont_per_bin(bin_contig_names, seed_contig_dict, cont_weight, marker_set_name='bacar',
                          cal_cont_strategy='checkm', rm_marker_minmax=False):
    """
    contigs_in_a_bin: contigs_list : the list contains the names of the contigs in a bin
    :return: # identified_marker_num, # uniq_marker_gene, comp, cont, comp -5*cont
    comp = #uniq_marker_gene_num/|marker gene set|
    cont = (identified_marker_num - uniq_marker_gene_num) /|marker gene set|
    cc_score = comp - 5*cont
    seed_contig_dict: key:marker name, value: list
    """
    identified_marker_num = 0
    uniq_marker_gene_num = 0

    if rm_marker_minmax:
        if marker_set_name == 'bacar':
            marker_set_size = 38
        elif marker_set_name == 'bacteria':
            marker_set_size = 103
        else:
            print("Please input the right marker_set_name.")
            sys.exit(1)
        remove_markers = ['gsa_f1k.fa.bacar.seed.0', 'gsa_f1k.fa.bacar.seed.39',
                          'gsa_f1k.fa.marker.seed.0', 'gsa_f1k.fa.marker.seed.1',
                          'gsa_f1k.fa.marker.seed.105', 'gsa_f1k.fa.marker.seed.106']
        for key in seed_contig_dict:
            if key not in remove_markers:
                identified_marker_num += len(list(set(seed_contig_dict[key]) & set(bin_contig_names)))
                if (len(list(set(seed_contig_dict[key]) & set(bin_contig_names))) > 0):
                    uniq_marker_gene_num += 1
    else:
        if marker_set_name == 'bacar':
            marker_set_size = 40
        elif marker_set_name == 'bacteria':
            marker_set_size = 107
        else:
            print("Please input the right marker_set_name.")
            sys.exit(1)

        for key in seed_contig_dict:
            identified_marker_num += len(list(set(seed_contig_dict[key]) & set(bin_contig_names)))
            if (len(list(set(seed_contig_dict[key]) & set(bin_contig_names))) > 0):
                uniq_marker_gene_num += 1

    comp = (uniq_marker_gene_num / marker_set_size) * 100

    if cal_cont_strategy == 'checkm':
        cont = 100 * ((identified_marker_num - uniq_marker_gene_num) / marker_set_size)
    elif cal_cont_strategy == 'amber':
        cont = 100 * ((identified_marker_num - uniq_marker_gene_num) / max(1, identified_marker_num))

    cc_score = comp - cont_weight * cont
    return identified_marker_num, uniq_marker_gene_num, comp, cont, cc_score


def select_comp_cont_per_bin(bin_contig_names, bacar_seed_contig_dict, bacteria_seed_contig_dict, cont_weight,
                             cal_cont_strategy='checkm', rm_marker_minmax=False):
    bacar_identified_marker_num, bacar_uniq_marker_gene_num, bacar_comp, bacar_cont, bacar_cc_score = cal_comp_cont_per_bin(
        bin_contig_names, bacar_seed_contig_dict,
        cont_weight, marker_set_name='bacar', cal_cont_strategy=cal_cont_strategy, rm_marker_minmax=rm_marker_minmax)
    bacteria_identified_marker_num, bacteria_uniq_marker_gene_num, bacteria_comp, bacteria_cont, bacteria_cc_score = cal_comp_cont_per_bin(
        bin_contig_names, bacteria_seed_contig_dict,
        cont_weight, marker_set_name='bacteria', cal_cont_strategy=cal_cont_strategy, rm_marker_minmax=rm_marker_minmax)

    if bacteria_cc_score > bacar_cc_score:
        domain = 'bacteria'
        return bacteria_identified_marker_num, bacteria_uniq_marker_gene_num, bacteria_comp, bacteria_cont, bacteria_cc_score, domain

    else:
        domain = 'bacar'
        return bacar_identified_marker_num, bacar_uniq_marker_gene_num, bacar_comp, bacar_cont, bacar_cc_score, domain


def gen_contig_list(fastafile):
    seq_ids = []
    with open(fastafile, 'r') as f:
        for line in f:
            if line.startswith(">"):
                if " " in line:
                    seq, others = line.split(' ', 1)
                    seq = line.rstrip("\n")
                    seq = seq.lstrip(">")
                    seq_ids.append(seq)
                else:
                    seq = line.rstrip("\n")
                    seq = seq.lstrip(">")
                    seq_ids.append(seq)
    return seq_ids


#############

############
def get_bin_extension(bin_dir):
    """Determine extension of bins."""
    exts = []
    for f in os.listdir(bin_dir):
        f_split = f.split('.')

        if len(f_split) > 1:
            ext = f_split[-1]
            if ext in ['fa', 'fna', 'fasta', 'bin']:
                exts.append(ext)
    if len(exts) == 0:
        return None, None

    ext, count = Counter(exts).most_common(1)[0]
    return ext, count


# changed from unitem
def get_bin_dirs(bin_dirs_file):
    bin_dirs = {}
    check_file_exists(bin_dirs_file)
    for line in open(bin_dirs_file):
        if line.strip():
            line_split = list(map(str.strip, line.split('\t')))
            if len(line_split) != 2:
                print("Skipping invalid line: %s" % line.strip())
                continue

            method_id = line_split[0]
            d = line_split[1]
            check_dir_exists(d)
            bin_ext, count = get_bin_extension(d)
            if not bin_ext:
                print('No bins identified for %s in %s.' % (method_id, d))
            else:
                bin_dirs[method_id] = (d, bin_ext)
                print("Processing %d genomes from %s with extension '%s'." % (count, method_id, bin_ext))
    return bin_dirs


# from unitem
# https://github.com/dparks1134/UniteM/blob/master/unitem/common.py
def read_bins(bin_dirs):
    """Read sequences in bins.
    bins : d[binning method][bin ID] -> contig_id
    Contigs for bins across all binning methods.
    contigs : d[cid] -> seq
    contigs_in_bins[seq_id][method_id]: d[contig_id][method_id]-> bin_id
    """

    bins = defaultdict(lambda: defaultdict(set))
    contigs = {}
    contigs_in_bins = defaultdict(lambda: {})
    for method_id, (bin_dir, bin_ext) in bin_dirs.items():
        for bf in os.listdir(bin_dir):
            if not bf.endswith(bin_ext):
                continue

            bin_id = bf[0:bf.rfind(bin_ext)]
            if bin_id[-1] == '.':
                bin_id = bin_id[0:-1]
            bf_path = os.path.join(bin_dir, bf)

            for seq_id, seq in seq_io.read_seq(bf_path):
                bins[method_id][bin_id].add(seq_id)
                contigs[seq_id] = seq
                contigs_in_bins[seq_id][method_id] = bin_id

            if len(bins[method_id][bin_id]) == 0:
                print('Bin %s from %s is empty.' % (bf, method_id))

    return bins, contigs, contigs_in_bins
