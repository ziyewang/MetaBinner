#!/usr/bin/env python

# -*- coding: utf-8 -*-
# @Author  : ZWang
# @FileName: component_binning.py (modified for Metabinner.py from earlier version)
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
from sklearn.cluster._kmeans import euclidean_distances, stable_cumsum, KMeans, check_random_state, row_norms, MiniBatchKMeans

from scipy.sparse import csc_matrix

logger = logging.getLogger('Metabinner v1.4.4')

logger.setLevel(logging.INFO)

# logging
formatter = logging.Formatter('%(asctime)s - %(message)s')

console_hdr = logging.StreamHandler()
console_hdr.setFormatter(formatter)

logger.addHandler(console_hdr)


def arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--contig_file', type=str, help=("The contigs file."))
    parser.add_argument('--coverage_profiles', type=str, help=(
        "The coverage profiles, containing a table where each row correspond to a contig, and each column correspond to a sequencing sample. All values are separated with tabs."))
    parser.add_argument('--composition_profiles', type=str, help=(
        "The composition profiles, containing a table where each row correspond to a contig, and each column correspond to the kmer composition of particular kmer. All values are separated with comma."))
    parser.add_argument('--output', type=str, help="The output file, storing the binning result.")
    parser.add_argument('--log', type=str, help="Specify where to store log file")
    parser.add_argument('--clusters', default=0, type=int,
                        help="Specify the number of clusters. If not specified, the cluster number is estimated by single-copy genes. If the specified number is smaller than the bin number estimated by MetaBinner, the cluster number will be determined as bin number estimated by MetaBinner.")
    parser.add_argument('--estimated_k', default=0, type=int,
                        help="Specify the number of estimated clusters by metabinner (only be used when users have obtained the bin number by metabinner). If specified, bin number estimation step will be skipped.")
    parser.add_argument('--threads', default=20, type=int,
                        help="the number of threads. default is 20.")
    parser.add_argument('--contig_length_threshold', default=1001, type=int,
                        help="The threshold of contig length for marker gene. default is 1001.")
    parser.add_argument('--dataset_scale', type=str, default="large", help=(
        "The scale of the dataset (for bin number identification),large or small, default is large. The parameter will affect the bin number estimation."))

    args = parser.parse_args()
    if not (args.contig_file and args.coverage_profiles and args.composition_profiles and args.output):
        parser.error(
            "Data is missing, add file(s) using --contig_file <contig_file> and/or --coverage_profiles <abund_profiles> and/or --composition_profiles <comp_profiles> and/or --output <out_file>")
        sys.exit(0)
    return args


def gen_X(com_file, cov_file):
    covHeader = pd.read_csv(cov_file, sep='\t', nrows=1)
    covMat = pd.read_csv(cov_file, sep='\t', usecols=range(1, covHeader.shape[1])).values
    namelist = pd.read_csv(cov_file, sep='\t', usecols=range(1)).values[:, 0]
    mapObj = dict(zip(namelist, range(len(namelist))))

    compositHeader = pd.read_csv(com_file, sep=',', nrows=1)
    shuffled_compositMat = pd.read_csv(com_file, sep=',', usecols=range(1, compositHeader.shape[1])).values
    shuffled_namelist = pd.read_csv(com_file, sep=',', usecols=range(1)).values[:, 0]

    covIdxArr = np.empty(len(mapObj), dtype=np.int)
    for contigIdx in range(len(shuffled_namelist)):
        if shuffled_namelist[contigIdx] in mapObj:
            covIdxArr[mapObj[shuffled_namelist[contigIdx]]] = contigIdx
    compositMat = shuffled_compositMat[covIdxArr]

    if covMat.shape[1] > 1:
        covMat = covMat + 1e-2
        covMat = covMat / covMat.sum(axis=0)[None, :]
        covMat = covMat / covMat.sum(axis=1)[:, None]
    else:
        covMat = covMat + 1e-2
        covMat = covMat / covMat.max(axis=0)[None, :]

    compositMat = compositMat + 1
    compositMat = compositMat / compositMat.sum(axis=1)[:, None]
    X_t = np.hstack((covMat, compositMat))  # del * 1e1
    return X_t, namelist, mapObj, covMat, compositMat


def gen_seed(contig_file, threads, contig_length_threshold, marker_name="marker", quarter="3quarter"):
    fragScanURL = 'run_FragGeneScan.pl'

    hmmExeURL = 'hmmsearch'
    markerExeURL = os.path.join(os.getcwd(), '../auxiliary', 'test_getmarker_' + quarter + '.pl')
    markerURL = os.path.join(os.getcwd(), '../auxiliary', marker_name + '.hmm')
    seedURL = contig_file + "." + marker_name + "." + quarter + "_lencutoff_" + str(contig_length_threshold) + ".seed"
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
                markerCmd = markerExeURL + " " + hmmResultURL + " " + contig_file + " " + str(
                    contig_length_threshold) + " " + seedURL
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
        if k < len(len_weight):
            if dataset_scale=='huge':
                kmeans = MiniBatchKMeans(n_clusters=k, random_state=7, n_init=30, init_size=min(len(len_weight),max(3*1024,3*k)))
            else:
                kmeans = KMeans(n_clusters=k, init='k-means++', random_state=7, n_init=30, n_jobs=threads)
            kmeans.fit(X_mat, sample_weight=len_weight)
            silVal = silhouette(X_mat, kmeans.cluster_centers_, kmeans.labels_, len_weight)
            logger.info("k:" + str(k) + "\tsilhouette:" + str(silVal) + "\telapsed time:" + str(time.time() - t))
            t = time.time()

            if silVal > bestSilVal:
                bestSilVal = silVal
                bestK = k
            else:
                break
        else:
            break
    candK = bestK + 2 * stepK
    bestSilVal_2nd = 0
    for k in range(candK, maxK, stepK):
        if k < len(len_weight):
            if dataset_scale=='huge':
                kmeans = MiniBatchKMeans(n_clusters=k, random_state=7, n_init=30, init_size=min(len(len_weight),max(3*1024,3*k)))
            else:
                kmeans = KMeans(n_clusters=k, init='k-means++', random_state=7, n_init=30, n_jobs=threads)
            kmeans.fit(X_mat, sample_weight=len_weight)
            silVal_2nd = silhouette(X_mat, kmeans.cluster_centers_, kmeans.labels_, len_weight)
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


def seed_kmeans_combo(seed_idx, output, X_mat, bin_number, prefix, length_weight, marker_name="marker1",
                      quarter="3quarter",dataset_scale="large", threads=-1):
    # run partial seed kmeans marker1_seed length weight
    logger.info("run partial seed kmeans " + marker_name + " seed length weight with:\t" + quarter + '_' + prefix)
    output_temp = os.path.dirname(
        output) + '/intermediate_result' + '/partial_seed_kmeans_' + marker_name + '_seed_length_weight_' + quarter + '_' + prefix + '_result.tsv'
    if not (os.path.exists(output_temp)):
        if dataset_scale == "huge":
            km = KMeans(n_clusters=bin_number, n_jobs=threads, random_state=7, n_init=30, algorithm="full",
                        init=functools.partial(partial_seed_init, seed_idx=seed_idx))
        else:
            km = KMeans(n_clusters=bin_number, n_jobs=threads, random_state=7, n_init=30,
                    init=functools.partial(partial_seed_init, seed_idx=seed_idx))
        km.fit(X_mat, sample_weight=length_weight)
        idx = km.labels_
        save_result(idx, output_temp, namelist)


def my_kmeans(X_mat, namelist, bin_number, bacar_marker_seed_num, length_weight, output, contig_length_threshold,
              prefix='X_t_notrans', quarter="3quarter", contig_file=None,dataset_scale='large',threads=-1):
    if bacar_marker_seed_num > 0:
        seed_bacar_marker_url = contig_file + ".bacar_marker" + "." + quarter + "_lencutoff_" + str(
            contig_length_threshold) + ".seed"
        seed_bacar_marker_idx = gen_seed_idx(seed_bacar_marker_url, contig_id_list=namelist)

    # run kmeans length weight
    logger.info("run kmeans length weight with:\t" + prefix)
    output_temp = os.path.dirname(output) + '/intermediate_result' + '/kmeans_length_weight_' + prefix + '_result.tsv'
    if not (os.path.exists(output_temp)):
        if dataset_scale == 'huge':
            km = MiniBatchKMeans(n_clusters=bin_number, random_state=7, n_init=30, init_size=min(len(length_weight),max(3*1024,3*bin_number)))
        else:
            km = KMeans(n_clusters=bin_number, init='k-means++', n_jobs=threads, n_init=30, random_state=7)
        km.fit(X_mat, sample_weight=length_weight)  # add log transform
        idx = km.labels_
        save_result(idx, output_temp, namelist)

    if bacar_marker_seed_num > 0:
        seed_kmeans_combo(seed_bacar_marker_idx, output, X_mat, bin_number, prefix, length_weight,
                          marker_name="bacar_marker", quarter=quarter,dataset_scale=dataset_scale, threads=threads)


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
    logger.info("Output file path:\t" + args.output)
    logger.info("Predefined Clusters:\t" + (str(args.clusters) if args.clusters > 0 else "Auto"))
    logger.info("The number of threads:\t" + str(args.threads))

    com_file = args.composition_profiles
    cov_file = args.coverage_profiles

    X_t, namelist, mapObj, X_cov, X_com = gen_X(com_file, cov_file)

    contigNum = X_t.shape[0]
    contig_file = args.contig_file
    output = args.output

    logger.info("The number of contigs:\t" + str(contigNum))

    clusters = args.clusters
    threads = args.threads
    contig_length_threshold = args.contig_length_threshold

    logger.info("gen bacar marker seed")
    bacar_marker_1quarter_seed_num = gen_seed(contig_file, threads, contig_length_threshold, marker_name="bacar_marker",
                                              quarter="1quarter")
    logger.info("bacar_marker_1quarter_seed_num:\t" + str(bacar_marker_1quarter_seed_num))
    bacar_marker_2quarter_seed_num = gen_seed(contig_file, threads, contig_length_threshold, marker_name="bacar_marker",
                                              quarter="2quarter")
    logger.info("bacar_marker_2quarter_seed_num:\t" + str(bacar_marker_2quarter_seed_num))
    bacar_marker_3quarter_seed_num = gen_seed(contig_file, threads, contig_length_threshold, marker_name="bacar_marker",
                                              quarter="3quarter")
    logger.info("bacar_marker_3quarter_seed_num:\t" + str(bacar_marker_3quarter_seed_num))

    logger.info("start calculate contig length")
    lengths = get_length(contig_file)
    length_weight = []
    for seq_id in namelist:
        length_weight.append(lengths[seq_id])

    dataset_scale = args.dataset_scale
    logger.info("Dataset scale:\t" + dataset_scale)

    # if dataset_scale == 'huge':
    #     X_t = np.log(X_t)
    #     X_t[np.abs(X_t)<1e-10]=0
    #     X_t = csc_matrix(X_t)
    #     X_cov = np.log(X_cov)
    #     X_cov[np.abs(X_cov)<1e-10]=0
    #     X_cov = csc_matrix(X_cov)
    #     X_com = np.log(X_com)
    #     X_com[np.abs(X_com)<1e-10]=0
    #     X_com = csc_matrix(X_com)
    # else:
    X_t = np.log(X_t)
    X_cov = np.log(X_cov)
    X_com = np.log(X_com)

    # set k0 using the large seed number from the two single-copy marker sets
    if args.estimated_k:
        bin_number = args.estimated_k
    else:
        # candK = max(marker1_3quarter_seed_num, bacar_marker_3quarter_seed_num) + 1
        candK = bacar_marker_3quarter_seed_num + 1
        logger.info("start estimate_bin_number")
        bin_number = estimate_bin_number(X_t, candK, dataset_scale=dataset_scale, len_weight=length_weight, threads=threads)
        logger.info("estimated_bin_number:\t" + str(bin_number))

    if args.clusters:
        # args.clusters is adopted only when it is greater than the bin_number
        bin_number = max(clusters, bin_number)

    intermediate_result_dir = os.path.dirname(output) + '/intermediate_result'
    if not (os.path.exists(intermediate_result_dir)):
        os.mkdir(intermediate_result_dir)
        # prefix shows the feature for binning
        # X_t



    my_kmeans(X_t, namelist, bin_number, bacar_marker_1quarter_seed_num, length_weight, output,
                  contig_length_threshold,
                  prefix='X_t_logtrans', quarter="1quarter", contig_file=contig_file,dataset_scale=dataset_scale, threads=threads)

    my_kmeans(X_t, namelist, bin_number, bacar_marker_2quarter_seed_num, length_weight, output,
                  contig_length_threshold,
                  prefix='X_t_logtrans', quarter="2quarter", contig_file=contig_file,dataset_scale=dataset_scale, threads=threads)

    my_kmeans(X_t, namelist, bin_number, bacar_marker_3quarter_seed_num, length_weight, output,
                  contig_length_threshold,
                  prefix='X_t_logtrans', quarter="3quarter", contig_file=contig_file, dataset_scale=dataset_scale, threads=threads)

    # X_com
    my_kmeans(X_com, namelist, bin_number, bacar_marker_1quarter_seed_num, length_weight, output,
                  contig_length_threshold,
                  prefix='X_com_logtrans', quarter="1quarter", contig_file=contig_file, dataset_scale=dataset_scale, threads=threads)

    my_kmeans(X_com, namelist, bin_number, bacar_marker_2quarter_seed_num, length_weight, output,
                  contig_length_threshold,
                  prefix='X_com_logtrans', quarter="2quarter", contig_file=contig_file, dataset_scale=dataset_scale, threads=threads)

    my_kmeans(X_com, namelist, bin_number, bacar_marker_3quarter_seed_num, length_weight, output,
                  contig_length_threshold,
                  prefix='X_com_logtrans', quarter="3quarter", contig_file=contig_file, dataset_scale=dataset_scale, threads=threads)

    # X_cov
    #if len(X_cov[0]) > 5:
    my_kmeans(X_cov, namelist, bin_number, bacar_marker_1quarter_seed_num,
                  length_weight, output, contig_length_threshold,
                  prefix='X_cov_logtrans', quarter="1quarter", contig_file=contig_file, dataset_scale=dataset_scale, threads=threads)

    my_kmeans(X_cov, namelist, bin_number, bacar_marker_2quarter_seed_num,
                  length_weight, output, contig_length_threshold,
                  prefix='X_cov_logtrans', quarter="2quarter", contig_file=contig_file, dataset_scale=dataset_scale, threads=threads)

    my_kmeans(X_cov, namelist, bin_number, bacar_marker_3quarter_seed_num,
                  length_weight, output, contig_length_threshold,
                  prefix='X_cov_logtrans', quarter="3quarter", contig_file=contig_file, dataset_scale=dataset_scale, threads=threads)
