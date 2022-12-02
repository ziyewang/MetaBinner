#!/usr/bin/env python 

# -*- coding: utf-8 -*-
# @Author  : ZWang
# @FileName: Metabinner.py
# scikit-learn == 0.22.1
# python 3.7

import numpy as np
import pandas as pd
import functools
import sys
import time
import gzip
import mimetypes
import os

import logging
import argparse

from Bio import SeqIO
import scipy.sparse as sp
from sklearn.cluster.k_means_ import euclidean_distances, stable_cumsum, KMeans, check_random_state, row_norms

from unitem_markers import Markers
from metabinner_util import get_bin_extension
from collections import defaultdict

from component_binning import gen_X

import biolib.seq_io as seq_io

logger = logging.getLogger('Metabinner post process for the component results')

logger.setLevel(logging.INFO)

# logging
formatter = logging.Formatter('%(asctime)s - %(message)s')

console_hdr = logging.StreamHandler()
console_hdr.setFormatter(formatter)

logger.addHandler(console_hdr)


# update for metabinner post process
def arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--contig_file', type=str, help=("The contigs file."))
    parser.add_argument('--coverage_profiles', type=str, help=(
        "The coverage profiles, containing a table where each row correspond to a contig, and each column correspond to a sample. All values are separated with tabs."))
    parser.add_argument('--composition_profiles', type=str, help=(
        "The composition profiles, containing a table where each row correspond to a contig, and each column correspond to the kmer composition of particular kmer. All values are separated with comma."))
    parser.add_argument('--ori_result_path', type=str, help="The path of the bins to be handled.")
    parser.add_argument('--bac_mg_table', type=str, help="The file of bac_mg_table.")
    parser.add_argument('--ar_mg_table', type=str, help="The file of ar_mg_table.")
    parser.add_argument('--log', type=str, help="Specify where to store log file")
    parser.add_argument('--threads', default=20, type=int,
                        help="the number of threads. default is 20.")
    parser.add_argument('--mincomp', default=70, type=int,
                        help="the mininum comp for post process. default is 70.")
    parser.add_argument('--mincont', default=50, type=int,
                        help="the mininum cont for post process. default is 50.")
    parser.add_argument('--dataset_scale', type=str, default="large", help=(
        "The scale of the dataset (for bin number identification),large or small, default is large"))

    args = parser.parse_args()
    if not (
            args.contig_file and args.coverage_profiles and args.composition_profiles and args.ori_result_path and args.bac_mg_table and args.ar_mg_table):
        parser.error(
            "Data is missing, add file(s) using --contig_file <contig_file> and/or --coverage_profiles <abund_profiles> and/or --composition_profiles <comp_profiles> and/or --ori_result_path <out_file> and/or --bac_mg_table and/or ar_mg_table")
        sys.exit(0)
    return args

#
# def gen_X(com_file, cov_file):
#     covHeader = pd.read_csv(cov_file, sep='\t', nrows=1)
#     covMat = pd.read_csv(cov_file, sep='\t', usecols=range(1, covHeader.shape[1])).values
#     namelist = pd.read_csv(cov_file, sep='\t', usecols=range(1)).values[:, 0]
#     mapObj = dict(zip(namelist, range(len(namelist))))
#
#     compositHeader = pd.read_csv(com_file, sep=',', nrows=1)
#     shuffled_compositMat = pd.read_csv(com_file, sep=',', usecols=range(1, compositHeader.shape[1])).values
#     shuffled_namelist = pd.read_csv(com_file, sep=',', usecols=range(1)).values[:, 0]
#
#     covIdxArr = np.empty(len(mapObj), dtype=np.int)
#     for contigIdx in range(len(shuffled_namelist)):
#         if shuffled_namelist[contigIdx] in mapObj:
#             covIdxArr[mapObj[shuffled_namelist[contigIdx]]] = contigIdx
#     compositMat = shuffled_compositMat[covIdxArr]
#
#     covMat = covMat + 1e-2
#     covMat = covMat / covMat.sum(axis=0)[None, :]
#     if covMat.shape[1] > 1:
#         covMat = covMat / covMat.sum(axis=1)[:, None]
#     compositMat = compositMat + 1
#     compositMat = compositMat / compositMat.sum(axis=1)[:, None]
#     X_t = np.hstack((covMat, compositMat))  # del * 1e1
#     return X_t, namelist, mapObj, covMat, compositMat

def gen_seed(contig_file, threads, marker_name="marker", quarter="3quarter"):
    fragScanURL = 'run_FragGeneScan.pl'

    hmmExeURL = 'hmmsearch'
    markerExeURL = os.path.join(os.getcwd(), 'auxiliary', 'test_getmarker_' + quarter + '.pl')
    markerURL = os.path.join(os.getcwd(), 'auxiliary', marker_name + '.hmm')
    seedURL = contig_file + "." + marker_name + "." + quarter + ".seed"
    fragResultURL = contig_file + ".frag.faa"
    hmmResultURL = contig_file + '.' + marker_name + ".hmmout"

    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL + " -genome=" + contig_file + " -out=" + contig_file + ".frag -complete=0 -train=complete -thread=" + str(
            threads) + " 1>" + contig_file + ".frag.out 2>" + contig_file + ".frag.err"
        logger.info("exec cmd: " + fragCmd)
        os.system(fragCmd)

    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL + " --domtblout " + hmmResultURL + " --cut_tc --cpu " + str(
                threads) + " " + markerURL + " " + fragResultURL + " 1>" + hmmResultURL + ".out 2>" + hmmResultURL + ".err"
            logger.info("exec cmd: " + hmmCmd)
            os.system(hmmCmd)

        if os.path.exists(hmmResultURL):
            if not (os.path.exists(seedURL)):
                markerCmd = markerExeURL + " " + hmmResultURL + " " + contig_file + " 1001 " + seedURL
                logger.info("exec cmd: " + markerCmd)
                os.system(markerCmd)

            if os.path.exists(seedURL):
                candK = file_len(seedURL)
            else:
                logger.info("markerCmd failed! Not exist: " + markerCmd)
                candK = 0
        else:
            logger.info("Hmmsearch failed! Not exist: " + hmmResultURL)
            sys.exit()
    else:
        logger.info("FragGeneScan failed! Not exist: " + fragResultURL)
        sys.exit()
    return candK

# estimate bin_number from candk
def estimate_bin_number(X_mat, candK, dataset_scale="large", len_weight=None,threads=-1):
    if dataset_scale == "small":
        candK = max(candK, 2)
        maxK = 4 * candK
        stepK = 2
    else:
        candK = max(candK, 2)
        maxK = 3 * candK
        stepK = 5
    bestK = candK
    bestSilVal = 0
    t = time.time()
    for k in range(candK, maxK, stepK):
        if k < len(X_mat):
            kmeans = KMeans(n_clusters=k, init='k-means++', random_state=7, n_init=30, n_jobs=threads)
            kmeans.fit(np.log(X_mat), sample_weight=len_weight)
            silVal = silhouette(np.log(X_mat), kmeans.cluster_centers_, kmeans.labels_, len_weight)
            logger.info("k:" + str(k) + "\tsilhouette:" + str(silVal) + "\telapsed time:" + str(time.time() - t))
            t = time.time()

            if silVal > bestSilVal:
                bestSilVal = silVal
                bestK = k
            else:
                break
        else:
            logger.info("bestk:" + str(bestK) + "\tsilVal:" + str(bestSilVal))
            return bestK
    candK = bestK + 2 * stepK
    bestSilVal_2nd = 0
    for k in range(candK, maxK, stepK):
        if k < len(X_mat):
            kmeans = KMeans(n_clusters=k, init='k-means++', random_state=7, n_init=30, n_jobs=threads)
            kmeans.fit(np.log(X_mat), sample_weight=len_weight)
            silVal_2nd = silhouette(np.log(X_mat), kmeans.cluster_centers_, kmeans.labels_, len_weight)
            logger.info("k:" + str(k) + "\tsilhouette:" + str(silVal_2nd) + "\telapsed time:" + str(time.time() - t))
            t = time.time()
            if silVal_2nd > bestSilVal_2nd:
                bestSilVal_2nd = silVal_2nd
                bestK = k
            else:
                break
        else:
            break
    if bestSilVal_2nd > bestSilVal:
        bestSilVal = bestSilVal_2nd
    else:
        bestK = candK - 2 * stepK
    logger.info("bestk:" + str(bestK) + "\tsilVal:" + str(bestSilVal))
    return bestK


def silhouette(X, W, label, len_weight):
    X_colsum = np.sum(X ** 2, axis=1)
    X_colsum = X_colsum.reshape(len(X_colsum), 1)
    W_colsum = np.sum(W ** 2, axis=1)
    W_colsum = W_colsum.reshape(len(W_colsum), 1)

    Dsquare = np.tile(X_colsum, (1, W.shape[0])) + np.tile(W_colsum.T, (X.shape[0], 1)) - 2 * X.dot(W.T)
    # avoid error caused by accuracy
    Dsquare[Dsquare < 0] = 0
    D = np.sqrt(Dsquare)
    aArr = D[np.arange(D.shape[0]), label]
    D[np.arange(D.shape[0]), label] = np.inf
    bArr = np.min(D, axis=1)
    tmp = (bArr - aArr) / np.maximum(aArr, bArr)
    return np.average(tmp, weights=len_weight)


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


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


def gen_seed_idx(seedURL, contig_id_list):
    seed_list = []
    with open(seedURL) as f:
        for line in f:
            if line.rstrip('\n') in contig_id_list:
                seed_list.append(line.rstrip('\n'))
    name_map = dict(zip(contig_id_list, range(len(contig_id_list))))
    seed_idx = [name_map[seed_name] for seed_name in seed_list]
    return seed_idx


def save_result(result, filepath, namelist):
    filedir, filename = os.path.split(filepath)
    if not filename:
        filename = "result.tsv"
    if not os.path.exists(filedir):
        os.makedirs(filedir)
    f = open(filepath, 'w')
    for contigIdx in range(len(result)):
        f.write(namelist[contigIdx] + "\t" + str(result[contigIdx].item(0)) + "\n")
    f.close()


# change from sklearn.cluster.kmeans
def partial_seed_init(X, n_clusters, random_state, seed_idx, n_local_trials=None):
    print('Using partial seed')

    random_state = check_random_state(random_state)
    x_squared_norms = row_norms(X, squared=True)

    n_samples, n_features = X.shape

    centers = np.empty((n_clusters, n_features), dtype=X.dtype)

    # Set the number of local seeding trials if none is given
    if n_local_trials is None:
        # This is what Arthur/Vassilvitskii tried, but did not report
        # specific results for other than mentioning in the conclusion
        # that it helped.
        n_local_trials = 2 + int(np.log(n_clusters))

    # Pick first center randomly

    center_id = seed_idx[0]

    if sp.issparse(X):
        centers[0] = X[center_id].toarray()
    else:
        centers[0] = X[center_id]

    # Initialize list of closest distances and calculate current potential
    closest_dist_sq = euclidean_distances(
        centers[0, np.newaxis], X, Y_norm_squared=x_squared_norms,
        squared=True)

    for c, center_id in enumerate(seed_idx[1:], 1):
        if sp.issparse(X):
            centers[c] = X[center_id].toarray()
        else:
            centers[c] = X[center_id]
        closest_dist_sq = np.minimum(closest_dist_sq,
                                     euclidean_distances(
                                         centers[c, np.newaxis], X, Y_norm_squared=x_squared_norms,
                                         squared=True))
    current_pot = closest_dist_sq.sum()

    # Pick the remaining n_clusters-1 points
    for c in range(len(seed_idx), n_clusters):
        # Choose center candidates by sampling with probability proportional
        # to the squared distance to the closest existing center
        rand_vals = random_state.random_sample(n_local_trials) * current_pot
        candidate_ids = np.searchsorted(stable_cumsum(closest_dist_sq),
                                        rand_vals)
        # XXX: numerical imprecision can result in a candidate_id out of range
        np.clip(candidate_ids, None, closest_dist_sq.size - 1,
                out=candidate_ids)

        # Compute distances to center candidates
        distance_to_candidates = euclidean_distances(
            X[candidate_ids], X, Y_norm_squared=x_squared_norms, squared=True)

        # Decide which candidate is the best
        best_candidate = None
        best_pot = None
        best_dist_sq = None
        for trial in range(n_local_trials):
            # Compute potential when including center candidate
            new_dist_sq = np.minimum(closest_dist_sq,
                                     distance_to_candidates[trial])
            new_pot = new_dist_sq.sum()

            # Store result if it is the best local trial so far
            if (best_candidate is None) or (new_pot < best_pot):
                best_candidate = candidate_ids[trial]
                best_pot = new_pot
                best_dist_sq = new_dist_sq

        # Permanently add best center candidate found in local tries
        if sp.issparse(X):
            centers[c] = X[best_candidate].toarray()
        else:
            centers[c] = X[best_candidate]
        current_pot = best_pot
        closest_dist_sq = best_dist_sq

    return centers


def save_result_refine(result, filepath, namelist, unclassified_contigs_id_number):
    filedir, filename = os.path.split(filepath)
    if not filename:
        filename = "result.tsv"
    if not os.path.exists(filedir):
        os.makedirs(filedir)
    f = open(filepath, 'w')
    for Idx in range(len(result)):
        f.write(namelist[unclassified_contigs_id_number[Idx]] + "\t" + str(result[Idx].item(0)) + "\n")
    f.close()


# convert 'bin' to 'fa'
def gen_bins(fastafile, resultfile, outputdir, prefix_str):
    # read fasta file
    logger.info("Processing file:\t{}".format(fastafile))
    sequences = {}
    if fastafile.endswith("gz"):
        with gzip.open(fastafile, 'r') as f:
            for line in f:
                line = str(line, encoding="utf-8")
                if line.startswith(">"):
                    if " " in line:
                        seq, others = line.split(' ', 1)
                        sequences[seq] = ""
                    else:
                        seq = line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    else:
        with open(fastafile, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    if " " in line:
                        seq, others = line.split(' ', 1)
                        sequences[seq] = ""
                    else:
                        seq = line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    logger.info("Reading Map:\t{}".format(resultfile))
    dic = {}
    with open(resultfile, "r") as f:
        for line in f:
            contig_name, cluster_name = line.strip().split('\t')  # change from split(',')
            try:
                dic[cluster_name].append(contig_name)
            except:
                dic[cluster_name] = []
                dic[cluster_name].append(contig_name)
    logger.info("Writing bins:\t{}".format(outputdir))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    bin_name = 0
    for _, cluster in dic.items():
        binfile = os.path.join(outputdir, "{}_{}.fa".format(prefix_str, bin_name))
        with open(binfile, "w") as f:
            for contig_name in cluster:
                contig_name = ">" + contig_name
                try:
                    sequence = sequences[contig_name]
                except:
                    bin_name += 1
                    continue
                f.write(contig_name + "\n")
                f.write(sequence + "\n")
                bin_name += 1


def read_bins_from_one_dir(bin_dir):
    bin_ext, count = get_bin_extension(bin_dir)
    bins = defaultdict(set)
    contigs = {}

    for bf in os.listdir(bin_dir):
        if not bf.endswith(bin_ext):
            continue
        bin_id = bf[0:bf.rfind(bin_ext)]
        if bin_id[-1] == '.':
            bin_id = bin_id[0:-1]
        bf_path = os.path.join(bin_dir, bf)

        for seq_id, seq in seq_io.read_seq(bf_path):
            bins[bin_id].add(seq_id)
            contigs[seq_id] = seq

    return bins, contigs


def split_hhbins(hhbin_contig_file, mapObj, X_t, length_weight, namelist, out_path, bin_id,threads=-1):
    hh_contigs_id = []
    for seq_record in SeqIO.parse(hhbin_contig_file, "fasta"):
        hh_contigs_id.append(seq_record.id)
    hh_contigs_id_number = [mapObj[x] for x in hh_contigs_id]
    X_t_hh_unclustered = X_t[hh_contigs_id_number]
    hh_weight = []
    for i in range(len(hh_contigs_id_number)):
        hh_weight.append(length_weight[hh_contigs_id_number[i]])

    seed_hh_num = gen_seed(hhbin_contig_file, threads, marker_name="bacar_marker")
    bin_number = estimate_bin_number(X_t_hh_unclustered, seed_hh_num, dataset_scale="small", len_weight=hh_weight,threads=threads)

    # seedurl may not exits??
    seedURL = hhbin_contig_file + ".bacar_marker.3_quarter.seed"
    # global seed_idx
    if os.path.exists(seedURL):
        seed_list = []
        with open(seedURL) as f:
            for line in f:
                if line.rstrip('\n') in namelist:
                    seed_list.append(line.rstrip('\n'))
        name_map = dict(zip(hh_contigs_id, range(len(hh_contigs_id))))
        seed_idx = [name_map[seed_name] for seed_name in seed_list]
        km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7,
                    init=functools.partial(partial_seed_init, seed_idx=seed_idx))
    else:
        km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7)

    km.fit(X_t_hh_unclustered, sample_weight=hh_weight)
    idx = km.labels_
    save_result_refine(idx, hhbin_contig_file + ".reclustered.tsv",
                       namelist, hh_contigs_id_number)
    gen_bins(hhbin_contig_file, hhbin_contig_file + ".reclustered.tsv",
             out_path, bin_id + "_reclustered")


if __name__ == '__main__':
    args = arguments()

    if args.log:
        handler = logging.FileHandler(args.log)
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    logger.info("Input arguments:")
    logger.info("Contig_file:\t" + args.contig_file)
    logger.info("Coverage_profiles:\t" + args.coverage_profiles)
    logger.info("Composition_profiles:\t" + args.composition_profiles)
    logger.info("The binning result file to be handled:\t" + args.ori_result_path)
    logger.info("The number of threads:\t" + str(args.threads))

    com_file = args.composition_profiles
    cov_file = args.coverage_profiles

    X_t, namelist, mapObj, X_cov, X_com = gen_X(com_file, cov_file)

    contigNum = X_t.shape[0]
    contig_file = args.contig_file

    logger.info("The number of contigs:\t" + str(contigNum))

    threads = args.threads

    markers = Markers()

    # bins, contigs = read_bins_from_one_dir(args.ori_result_path)
    path = args.ori_result_path
    bins, contigs = read_bins_from_one_dir(path)

    contig_lens = {cid: len(contigs[cid]) for cid in contigs}

    length_weight = []
    for seq_id in namelist:
        length_weight.append(contig_lens[seq_id])

    gene_tables = markers.marker_gene_tables(args.bac_mg_table, args.ar_mg_table)

    bin_dir = path
    out_path = bin_dir + '_post_process_mincomp_' + str(args.mincomp) + '_mincont_' + str(args.mincont) + '_bins/'
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    bin_ext, count = get_bin_extension(bin_dir)

    q = []
    for bf in os.listdir(bin_dir):
        if not bf.endswith(bin_ext):
            continue
        bin_id = bf[0:bf.rfind(bin_ext)]
        if bin_id[-1] == '.':
            bin_id = bin_id[0:-1]
        bf_path = os.path.join(bin_dir, bf)

        # for seq_id, seq in seq_io.read_seq(bf_path):
        #     bins[bin_id].add(seq_id)
        #     contigs[seq_id] = seq

        domain, comp, cont = markers.bin_quality(bins[bin_id])
        if comp >= float(args.mincomp) and cont >= float(args.mincont) and len(bins[bin_id]) >= 3:
            split_hhbins(bf_path, mapObj, X_t, length_weight, namelist, out_path, bin_id, threads=threads)
        else:
            temp_bin_file = os.path.join(out_path + bin_id + '.fa')
            fout_bin = open(temp_bin_file, 'w')
            for seq_id in bins[bin_id]:
                fout_bin.write('>%s\n' % seq_id)
                fout_bin.write(contigs[seq_id] + '\n')

            fout_bin.close()

        q.append((domain, comp, cont))
