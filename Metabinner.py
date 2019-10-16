# -*- coding: utf-8 -*-
# @Author  : ZWang
# @FileName: Metabinner.py
# scikit-learn == 0.20.4

import argparse
import csv
import logging
import math
import mimetypes
import os
import re
import shutil
import subprocess
import sys
import time
import gzip
import functools
from argparse import RawTextHelpFormatter

import numpy as np
import pandas as pd
import util
from Bio import SeqIO
# from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
import sklearn.cluster.k_means_ as kmeans
import scipy.sparse as sp
from sklearn.cluster.k_means_ import euclidean_distances, stable_cumsum, KMeans, check_random_state, row_norms
from sklearn.metrics import pairwise_distances

from multiprocessing import Pool
# from concurrent.futures import ProcessPoolExecutor
# pool = Pool()
from functools import partial

logger = logging.getLogger('Metabinner')

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
        "The coverage profiles, containing a table where each row correspond to a contig, and each column correspond to a sample. All values are separated with tabs."))
    parser.add_argument('--composition_profiles', type=str, help=(
        "The composition profiles, containing a table where each row correspond to a contig, and each column correspond to the kmer composition of particular kmer. All values are separated with comma."))
    parser.add_argument('--output', type=str, help="The output file, storing the binning result.")
    parser.add_argument('--log', type=str, help="Specify where to store log file")
    parser.add_argument('--clusters', default=0, type=int,
                        help="Specify the number of clusters. If not specified, the cluster number is estimated by single-copy genes.")
    """parser.add_argument('--weight_length', action="store_true",
                        help="use contig length as the weight of the contigs")  ##?
    parser.add_argument('--seed', action='store_true',
                        help="use part of the contigs with single-marker gene as the initial of the kmeans")"""
    parser.add_argument('--use_hmm', action="store_true", help="use hmm profile as another feature")
    parser.add_argument('--hmm_icm_path', type=str,
                        help="The path that contains trained imm model (please end with '/').")
    parser.add_argument('--hmm_file', type=str, help=("The hmm profiles."))
    parser.add_argument('--pacbio_read_profiles', type=str, help=(
        "The pacbio_read coverage profiles, containing a table where each row correspond to a contig, and each column correspond to a sample. All values are separated with tabs."))
    parser.add_argument('--binscore', default=0.3, type=float,
                        help="Specify the score threshold for das_tool.")
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

    covMat = covMat + 1e-2
    covMat = covMat / covMat.sum(axis=0)[None, :]
    if covMat.shape[1] >= 10:
        covMat = covMat / covMat.sum(axis=1)[:, None]
    compositMat = compositMat + 1
    compositMat = compositMat / compositMat.sum(axis=1)[:, None]
    X_t = np.hstack((covMat, compositMat)) * 1e1
    return X_t, namelist, mapObj, covMat, compositMat


def gen_X_cov_pb(cov_pb_file, mapObj):
    cov_pb_Header = pd.read_csv(cov_pb_file, sep='\t', nrows=1)
    shuffled_cov_pb_Mat = pd.read_csv(cov_pb_file, sep='\t', usecols=range(1, cov_pb_Header.shape[1])).values
    shuffled_namelist = pd.read_csv(cov_pb_file, sep='\t', usecols=range(1)).values[:, 0]

    covIdxArr = []
    for contigIdx in range(len(shuffled_namelist)):
        assert (shuffled_namelist[contigIdx] in mapObj)
        covIdxArr.append(mapObj[shuffled_namelist[contigIdx]])
    cov_pb_Mat = shuffled_cov_pb_Mat.copy()
    cov_pb_Mat[covIdxArr, :] = shuffled_cov_pb_Mat

    cov_pb_Mat = cov_pb_Mat + 1e-2
    cov_pb_Mat = cov_pb_Mat / cov_pb_Mat.sum(axis=0)[None, :]
    if cov_pb_Mat.shape[1] >= 10:
        cov_pb_Mat = cov_pb_Mat / cov_pb_Mat.sum(axis=1)[:, None]

    return cov_pb_Mat


def gen_X_hmm(hmm_file, cov_file):
    print("hmm_profiles:\t" + hmm_file)
    print("abundance_profiles:\t" + cov_file)
    namelist = pd.read_csv(cov_file, sep='\t', usecols=range(1)).values[:, 0]
    mapObj = dict(zip(namelist, range(len(namelist))))
    hmmHeader = pd.read_csv(hmm_file, sep='\t', nrows=1)
    shuffled_hmmMat = pd.read_csv(hmm_file, sep='\t', usecols=range(1, hmmHeader.shape[1])).values
    shuffled_namelist = pd.read_csv(hmm_file, sep='\t', usecols=range(1)).values[:, 0]
    covIdxArr = np.empty(len(mapObj), dtype=np.int)

    for contigIdx in range(len(shuffled_namelist)):  # 就这一处用到shuffled_namelist
        if shuffled_namelist[contigIdx] in mapObj:
            covIdxArr[mapObj[shuffled_namelist[contigIdx]]] = contigIdx

    hmmMat = shuffled_hmmMat[covIdxArr]

    return hmmMat


# changed from scimm
def score_reads(path, files, readsf, par, output):  # 需要改路径
    cmds = []
    hmm_score_dir = os.path.dirname(output) + '/hmm_score'
    simple_scoreURL = os.path.join(os.getcwd(), 'simple-score')
    os.system("chmod 777 " + simple_scoreURL)
    os.mkdir(hmm_score_dir)
    for file in files:
        """if not os.path.isdir(path + file):
            cmds.append(
                '/home/wzy/binning/scimm/bin/simple-score -N %shmm/%s.icm < %s > %sscore/%s.scores.tmp ' % (
                    path, file, readsf, path, file))"""
        if not os.path.isdir(path + file):
            cmds.append(
                '%s -N %s%s < %s > %s/%s.scores.tmp ' % (simple_scoreURL,
                                                         path, file, readsf, hmm_score_dir, file))

    util.exec_par(cmds, par)


def get_read_probs(path, files, output):
    hmm_score_dir = os.path.dirname(output) + '/hmm_score/'
    k = len(os.listdir(hmm_score_dir))
    read_likes = {}
    c = 0
    for file in files:
        if not os.path.isdir(path + file):
            for line in open('%s%s.scores.tmp' % (hmm_score_dir, file)):
                (r, s) = line.split('\t')
                r = r.strip()
                # 增加能去掉字符串k141_12130 flag=1 multi=29.9095 len=7194后面一�?
                r = r.split()[0]
                if not r in read_likes:
                    read_likes[r] = [0] * k
                read_likes[r][c] = float(s)
            c = c + 1
    read_probs = {}
    likelihood = 0.0
    for r in read_likes:
        # combine mate likelihoods and priors
        r1 = read_likes[r]
        # read_scores = [r1[x] + math.log(priors[x]) for x in range(k)]
        read_scores = [r1[x] for x in range(k)]
        # determine probabilities of assignments
        sum_score = read_scores[0]
        for i in range(1, k):
            sum_score = log_add(sum_score, read_scores[i])
        read_probs[r] = []
        for i in range(k):
            sc = math.exp(read_scores[i] - sum_score)
            if sc > 1e-8:
                read_probs[r].append(sc)
            else:
                read_probs[r].append(0)
        # update likelihood, accounting for mates being assigned twice
        likelihood += max(read_scores)
    return likelihood, read_probs


def log_add(l_i, l_j):
    if l_i > l_j:
        return l_i + math.log(1 + math.exp((l_j - l_i)))
    else:
        return l_j + math.log(1 + math.exp((l_i - l_j)))


def gen_bestk(contig_file, X_mat, bestK=0):  # 改成无论是否固定k都要跑生成seed，去掉没有生成的话从5开始随机的情况
    fragScanURL = os.path.join(os.getcwd(), 'auxiliary', 'FragGeneScan1.19', 'run_FragGeneScan.pl')
    os.system("chmod +777 " + fragScanURL)
    hmmExeURL = os.path.join(os.getcwd(), 'auxiliary', 'hmmer-3.1b1', 'bin', 'hmmsearch')
    os.system("chmod 777 " + hmmExeURL)
    markerExeURL = os.path.join(os.getcwd(), 'auxiliary', 'test_getmarker.pl')
    os.system("chmod 777 " + markerExeURL)
    markerURL = os.path.join(os.getcwd(), 'auxiliary', 'marker.hmm')
    seedURL = contig_file + ".seed"
    fragResultURL = contig_file + ".frag.faa"
    hmmResultURL = contig_file + ".hmmout"

    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL + " -genome=" + contig_file + " -out=" + contig_file + ".frag -complete=0 -train=complete -thread=48 1>" + contig_file + ".frag.out 2>" + contig_file + ".frag.err"
        logger.info("exec cmd: " + fragCmd)
        os.system(fragCmd)

    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL + " --domtblout " + hmmResultURL + " --cut_tc --cpu 48 " + markerURL + " " + fragResultURL + " 1>" + hmmResultURL + ".out 2>" + hmmResultURL + ".err"
            logger.info("exec cmd: " + hmmCmd)
            os.system(hmmCmd)

        if os.path.exists(hmmResultURL):
            if not (os.path.exists(seedURL)):
                markerCmd = markerExeURL + " " + hmmResultURL + " " + contig_file + " 1000 " + seedURL
                logger.info("exec cmd: " + markerCmd)
                os.system(markerCmd)

            if os.path.exists(seedURL):
                candK = file_len(seedURL)
                maxK = 2 * candK
                stepK = 2
            else:
                logger.info("seed not exist, k start from 2 ")
                # sys.exit()
                candK = 2
                maxK = 20
                stepK = 2
        else:
            logger.info("Hmmsearch failed! Not exist: " + hmmResultURL)
            sys.exit()
    else:
        logger.info("FragGeneScan failed! Not exist: " + fragResultURL)
        sys.exit()

    if bestK == 0:
        bestK = candK
        bestSilVal = 0
        t = time.time()
        for k in range(candK, maxK, stepK):
            kmeans = KMeans(n_clusters=k, init='k-means++', random_state=9, n_jobs=-1)
            kmeans.fit(X_mat)
            silVal = silhouette(X_mat, kmeans.cluster_centers_, kmeans.labels_)
            logger.info("k:" + str(k) + "\tsilhouette:" + str(silVal) + "\telapsed time:" + str(time.time() - t))
            t = time.time()

            if silVal > bestSilVal:
                bestSilVal = silVal
                bestK = k
            else:
                break

        candK = bestK + 4
        bestSilVal_2nd = 0
        for k in range(candK, maxK, stepK):
            kmeans = KMeans(n_clusters=k, init='k-means++', random_state=9, n_jobs=-1)
            kmeans.fit(X_mat)
            silVal_2nd = silhouette(X_mat, kmeans.cluster_centers_, kmeans.labels_)
            logger.info("k:" + str(k) + "\tsilhouette:" + str(silVal_2nd) + "\telapsed time:" + str(time.time() - t))
            t = time.time()
            if silVal_2nd > bestSilVal_2nd:
                bestSilVal_2nd = silVal_2nd
                bestK = k
            else:
                break
        if bestSilVal_2nd > bestSilVal:
            bestSilVal = bestSilVal_2nd
        else:
            bestK = candK - 4
        logger.info("bestk:" + str(bestK) + "\tsilVal:" + str(bestSilVal))

    else:
        logger.info("Use the pre-specified cluster number! k=" + str(bestK))

    return bestK


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


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def silhouette(X, W, label):
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
    return np.mean(tmp)


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
        binfile = os.path.join(outputdir, "{}_{}.bin".format(prefix_str, bin_name))
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


# change from binsanity
def checkm_analysis(file_, suffix_str, output):
    file_object = open(file_, 'r')
    lines = file_object.readlines()
    file_deal = open(file_ + '_deal.txt', 'w')
    try:
        for line in lines:
            if line.startswith(' '):
                print(line)
                file_deal.writelines(line)
    finally:
        file_object.close()
    file_deal.close()  # 要不会影响读�?

    checkm = list(csv.reader(open(file_ + '_deal.txt', 'r')))
    new = []
    for list_ in checkm:
        x = re.sub(' +', ' ', str(re.split(r'\t+', list_[0].rstrip('\t'))))
        new.append(x)

    del new[0]

    checkm_info_list = [list_.strip("['']") for list_ in new]
    checkm_info_list = [x.split() for x in checkm_info_list]

    good_bins = []
    High_completion_high_contamination = []
    # low_completion=[]
    others = []

    for list_ in checkm_info_list:
        if ((float(list_[12]) > 70 and (float(list_[13]) < 15)) or (float(list_[12]) - 5 * float(list_[13])) > 50):
            good_bins.append(list_[0])
        elif (float(list_[12]) > 70 and (float(list_[13]) > 50)):
            High_completion_high_contamination.append(list_[0])
        else:
            others.append(list_[0])

    if os.path.isdir(os.path.dirname(output) + "/good_bins") is False:
        os.makedirs(os.path.dirname(output) + "/good_bins")
    if os.path.isdir(os.path.dirname(output) + "/High_completion_high_contamination") is False:
        os.makedirs(os.path.dirname(output) + "/High_completion_high_contamination")
    if os.path.isdir(os.path.dirname(output) + "/others") is False:
        os.makedirs(os.path.dirname(output) + "/others")

    for name in good_bins:
        shutil.move((os.path.dirname(output) + '/' + str(name) + suffix_str),
                    os.path.dirname(output) + "/good_bins")
    for name in High_completion_high_contamination:
        shutil.move((os.path.dirname(output) + '/' + str(name) + suffix_str),
                    os.path.dirname(output) + "/High_completion_high_contamination")
    for name in others:
        shutil.move((os.path.dirname(output) + '/' + str(name) + suffix_str), os.path.dirname(output) + "/others")


# change from sklearn.cluster.kmeans
def partial_seed_init(X, n_clusters, random_state, seed_idx, n_local_trials=None):
    logger.info('Partial Seed Initialization')

    random_state = check_random_state(random_state)
    x_squared_norms = row_norms(X, squared=True)

    n_samples, n_features = X.shape
    centers = np.empty((n_clusters, n_features), dtype=X.dtype)
    assert x_squared_norms is not None, 'x_squared_norms None in _k_init'
    # Set the number of local seeding trials if none is given
    if n_local_trials is None:
        # This is what Arthur/Vassilvitskii tried, but did not report
        # specific results for other than mentioning in the conclusion
        # that it helped.
        n_local_trials = 2 + int(np.log(n_clusters))

    # global seed_idx
    # Pick first center randomly  ###修改�?
    if len(seed_idx) != n_samples:
        while True:
            center_id = random_state.randint(n_samples)
            if center_id not in seed_idx:
                break

        if sp.issparse(X):
            centers[0] = X[center_id].toarray()
        else:
            centers[0] = X[center_id]

        # Initialize list of closest distances and calculate current potential
        closest_dist_sq = euclidean_distances(
            centers[0, np.newaxis], X, Y_norm_squared=x_squared_norms,
            squared=True)

        current_pot = closest_dist_sq.sum()

        # Pick the remaining n_clusters-1 points

        for c in range(1, n_clusters - len(seed_idx)):  # changed
            # Choose center candidates by sampling with probability proportional

            # to the squared distance to the closest existing center

            rand_vals = random_state.random_sample(n_local_trials) * current_pot
            candidate_ids = np.searchsorted(stable_cumsum(closest_dist_sq),

                                            rand_vals)

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
                if ((best_candidate is None) or (new_pot < best_pot)) and (
                        (candidate_ids[trial]) not in seed_idx):  # changed
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
    if sp.issparse(X):
        centers[(n_clusters - len(seed_idx)):n_clusters] = X[seed_idx].toarray()
    else:
        centers[(n_clusters - len(seed_idx)):n_clusters] = X[seed_idx]
    # print(seed_idx)
    return centers


def recluster_other_contigs(not_clustered_path, X_t, namelist, mapObj, length_weight):
    files = os.listdir(not_clustered_path)
    other_contig_file = not_clustered_path + '/init_unclustered_contigs.fa'
    ofile = open(other_contig_file, 'w')
    # 遍历读取所有文件，并写入到输出文件
    for fr in files:
        if fr != 'init_unclustered_contigs.fa':
            for txt in open(not_clustered_path + '/' + fr, 'r'):
                ofile.write(txt)
    ofile.close()

    unclassified_contigs_id = []
    for seq_record in SeqIO.parse(other_contig_file, "fasta"):
        unclassified_contigs_id.append(seq_record.id)
    unclassified_contigs_id_number = [mapObj[x] for x in unclassified_contigs_id]
    X_t_unclustered = X_t[unclassified_contigs_id_number]

    bin_number = gen_bestk(other_contig_file, X_t_unclustered, 0)
    logger.info("bin_number for other contigs: %d", bin_number)

    unclassified_contigs_weight = []
    for i in range(len(unclassified_contigs_id_number)):
        unclassified_contigs_weight.append(length_weight[unclassified_contigs_id_number[i]])

    seedURL = other_contig_file + ".seed"
    # global seed_idx
    if os.path.exists(seedURL):
        seed_list = []
        with open(seedURL) as f:
            for line in f:
                seed_list.append(line.rstrip('\n'))
        name_map = dict(zip(unclassified_contigs_id, range(len(unclassified_contigs_id))))
        seed_idx = [name_map[seed_name] for seed_name in seed_list]
        km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7,
                    init=functools.partial(partial_seed_init, seed_idx=seed_idx))
    else:
        km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7)

    logger.info("Start bin the other bins.")
    # 之后可以加GC，像binsanity那样
    km.fit(X_t_unclustered, sample_weight=unclassified_contigs_weight)
    idx = km.labels_
    not_clustered_path_output = not_clustered_path + 'reclustered_result.tsv'
    save_result_refine(idx, not_clustered_path_output, namelist, unclassified_contigs_id_number)
    gen_bins(other_contig_file, not_clustered_path_output, os.path.dirname(not_clustered_path) + '/good_bins',
             "reclustered")


def gen_seed_number(contig_file):
    fragScanURL = os.path.join(os.getcwd(), 'auxiliary', 'FragGeneScan1.19', 'run_FragGeneScan.pl')
    os.system("chmod 777 " + fragScanURL)
    hmmExeURL = os.path.join(os.getcwd(), 'auxiliary', 'hmmer-3.1b1', 'bin', 'hmmsearch')
    os.system("chmod 777 " + hmmExeURL)
    markerExeURL = os.path.join(os.getcwd(), 'auxiliary', 'test_getmarker.pl')
    os.system("chmod 777 " + markerExeURL)
    markerURL = os.path.join(os.getcwd(), 'auxiliary', 'marker.hmm')
    seedURL = contig_file + ".seed"
    fragResultURL = contig_file + ".frag.faa"
    hmmResultURL = contig_file + ".hmmout"

    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL + " -genome=" + contig_file + " -out=" + contig_file + ".frag -complete=0 -train=complete -thread=10 1>" + contig_file + ".frag.out 2>" + contig_file + ".frag.err"
        logger.info("exec cmd: " + fragCmd)
        os.system(fragCmd)

    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL + " --domtblout " + hmmResultURL + " --cut_tc --cpu 10 " + markerURL + " " + fragResultURL + " 1>" + hmmResultURL + ".out 2>" + hmmResultURL + ".err"
            logger.info("exec cmd: " + hmmCmd)
            os.system(hmmCmd)

        if os.path.exists(hmmResultURL):
            if not (os.path.exists(seedURL)):
                markerCmd = markerExeURL + " " + hmmResultURL + " " + contig_file + " 1000 " + seedURL
                logger.info("exec cmd: " + markerCmd)
                os.system(markerCmd)

            if os.path.exists(seedURL):
                candK = file_len(seedURL)
            else:
                logger.info("markerCmd failed! Not exist: " + markerCmd)
                candK = 2
        else:
            logger.info("Hmmsearch failed! Not exist: " + hmmResultURL)
            sys.exit()
    else:
        logger.info("FragGeneScan failed! Not exist: " + fragResultURL)
        sys.exit()
    return candK


def gen_seed_number_bacar_marker(contig_file):
    fragScanURL = os.path.join(os.getcwd(), 'auxiliary', 'FragGeneScan1.19', 'run_FragGeneScan.pl')
    os.system("chmod 777 " + fragScanURL)
    hmmExeURL = os.path.join(os.getcwd(), 'auxiliary', 'hmmer-3.1b1', 'bin', 'hmmsearch')
    os.system("chmod 777 " + hmmExeURL)
    markerExeURL = os.path.join(os.getcwd(), 'auxiliary', 'test_getmarker.pl')
    os.system("chmod 777 " + markerExeURL)
    markerURL = os.path.join(os.getcwd(), 'auxiliary', 'bacar_marker.hmm')
    seedURL = contig_file + ".bacar_marker.seed"
    fragResultURL = contig_file + ".frag.faa"
    hmmResultURL = contig_file + ".bacar_marker.hmmout"

    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL + " -genome=" + contig_file + " -out=" + contig_file + ".frag -complete=0 -train=complete -thread=20 1>" + contig_file + ".frag.out 2>" + contig_file + ".frag.err"
        logger.info("exec cmd: " + fragCmd)
        os.system(fragCmd)

    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL + " --domtblout " + hmmResultURL + " --cut_tc --cpu 20 " + markerURL + " " + fragResultURL + " 1>" + hmmResultURL + ".out 2>" + hmmResultURL + ".err"
            logger.info("exec cmd: " + hmmCmd)
            os.system(hmmCmd)

        if os.path.exists(hmmResultURL):
            if not (os.path.exists(seedURL)):
                markerCmd = markerExeURL + " " + hmmResultURL + " " + contig_file + " 1000 " + seedURL
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

    if candK == 0:
        bestK = 0
    else:
        bestK = candK
        maxK = 2 * candK
        stepK = 2
        bestSilVal = 0
        t = time.time()
        for k in range(candK, maxK, stepK):
            kmeans = KMeans(n_clusters=k, init='k-means++', random_state=9, n_jobs=-1)
            kmeans.fit(X_t)
            silVal = silhouette(X_t, kmeans.cluster_centers_, kmeans.labels_)
            logger.info("k:" + str(k) + "\tsilhouette:" + str(silVal) + "\telapsed time:" + str(time.time() - t))
            t = time.time()

            if silVal > bestSilVal:
                bestSilVal = silVal
                bestK = k
            else:
                break

        candK = bestK + 4
        bestSilVal_2nd = 0
        for k in range(candK, maxK, stepK):
            kmeans = KMeans(n_clusters=k, init='k-means++', random_state=9, n_jobs=-1)
            kmeans.fit(X_t)
            silVal_2nd = silhouette(X_t, kmeans.cluster_centers_, kmeans.labels_)
            logger.info("k:" + str(k) + "\tsilhouette:" + str(silVal_2nd) + "\telapsed time:" + str(time.time() - t))
            t = time.time()
            if silVal_2nd > bestSilVal_2nd:
                bestSilVal_2nd = silVal_2nd
                bestK = k
            else:
                break
        if bestSilVal_2nd > bestSilVal:
            bestSilVal = bestSilVal_2nd
        else:
            bestK = candK - 4
        logger.info("bestk:" + str(bestK) + "\tsilVal:" + str(bestSilVal))

    return bestK


def recluster_hh_bins(high_com_p_high_cont_path, mapObj, X_t, length_weight, namelist):
    hh_files = os.listdir(high_com_p_high_cont_path)
    for file_ in hh_files:
        if file_.endswith('.bin'):
            hh_contigs_id = []
            hh_contig_file = high_com_p_high_cont_path + '/' + file_
            for seq_record in SeqIO.parse(hh_contig_file, "fasta"):
                hh_contigs_id.append(seq_record.id)
            hh_contigs_id_number = [mapObj[x] for x in hh_contigs_id]
            X_t_hh_unclustered = X_t[hh_contigs_id_number]
            bin_number = gen_bestk(hh_contig_file, X_t_hh_unclustered, 0)

            hh_weight = []
            for i in range(len(hh_contigs_id_number)):
                hh_weight.append(length_weight[hh_contigs_id_number[i]])
            # seedurl不一定存�?
            seedURL = hh_contig_file + ".seed"
            # global seed_idx
            if os.path.exists(seedURL):
                seed_list = []
                with open(seedURL) as f:
                    for line in f:
                        seed_list.append(line.rstrip('\n'))
                name_map = dict(zip(hh_contigs_id, range(len(hh_contigs_id))))
                seed_idx = [name_map[seed_name] for seed_name in seed_list]
                km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7,
                            init=functools.partial(partial_seed_init, seed_idx=seed_idx))
            else:
                km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7)

            km.fit(X_t_hh_unclustered, sample_weight=hh_weight)
            idx = km.labels_
            save_result_refine(idx, hh_contig_file + ".reclustered.tsv",
                               namelist, hh_contigs_id_number)
            gen_bins(hh_contig_file, hh_contig_file + ".reclustered.tsv",
                     os.path.dirname(high_com_p_high_cont_path) + '/good_bins', file_ + "_reclustered")


def read_fasta_file(fasta_file):
    with open(fasta_file, 'r') as read_handler:
        for line in read_handler:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                yield line[1:]


def convert(paths, output_file):
    files = os.listdir(paths)
    fasta_files = []
    for file in files:
        if file.endswith(('.fasta', '.fa', '.fna', '.bin')):
            fasta_files.append(file)
    with open(output_file, 'w') as write_handler:
        for bin_id, fasta_file in enumerate(fasta_files):
            for sequence_id in read_fasta_file(paths + '/' + fasta_file):
                write_handler.write("%s\t%s\n" % (sequence_id, bin_id))  # change from ","


def calculate_eps_for_l1(goodbin_path, mapObj, X_t, length_weight, namelist):
    goodbin_files = os.listdir(goodbin_path)
    score = []
    for file_ in goodbin_files:
        if file_.endswith('.bin'):
            goodbin_contigs_id = []
            goodbin_contig_file = goodbin_path + '/' + file_
            for seq_record in SeqIO.parse(goodbin_contig_file, "fasta"):
                goodbin_contigs_id.append(seq_record.id)
            goodbin_contigs_id_number = [mapObj[x] for x in goodbin_contigs_id]
            X_t_goodbin = X_t[goodbin_contigs_id_number]

            l1_distance_score = pairwise_distances(X_t_goodbin, metric='l1', n_jobs=-1)
            print(len(X_t_goodbin))
            if len(X_t_goodbin) >= 20:  # 这个20的数是根据minsamples的信息估算的
                l1_distance_score = np.sort(l1_distance_score, axis=-1)[:, 19]
                l1_distance_score = np.mean(l1_distance_score)
                # l1_distance_score=np.max(l1_distance_score)#这个想法不对，eps是邻�?
                score.append(l1_distance_score)
    return score


def gen_remained_fasta_file(fastafile, remained_contig_id, outputdir, prefix_str):
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
    dic = {}
    cluster_name = 'remained'
    for contig_name in remained_contig_id:
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
        binfile = os.path.join(outputdir, "{}_{}.bin".format(prefix_str, bin_name))
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


if __name__ == '__main__':
    args = arguments()
    # kmeans._k_init = _k_init

    if args.log:
        handler = logging.FileHandler(args.log)
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    logger.info("Input arguments:")
    logger.info("contig_file:\t" + args.contig_file)
    logger.info("coverage_profiles:\t" + args.coverage_profiles)
    logger.info("composition_profiles:\t" + args.composition_profiles)
    logger.info("output path:\t" + args.output)
    logger.info("clusters:\t" + (str(args.clusters) if args.clusters > 0 else "Auto"))

    com_file = args.composition_profiles
    cov_file = args.coverage_profiles

    X_t, namelist, mapObj, X_cov_sr, X_com = gen_X(com_file, cov_file)
    contigNum = X_t.shape[0]
    contig_file = args.contig_file

    if args.pacbio_read_profiles:
        X_cov_pb = gen_X_cov_pb(args.pacbio_read_profiles, mapObj)
        X_cov = np.hstack((X_cov_sr, X_cov_pb))
        X_t = np.hstack((X_t, X_cov_pb * 1e1))
    else:
        X_cov = X_cov_sr

    if args.use_hmm:
        if not args.hmm_file:
            path = args.hmm_icm_path  # '/mnt/data3/wzy/ncbi_genomes_metawatt_process/filter_small_genome/'  # 需要改路径
            files = os.listdir(path)
            par = 40
            # train_imm(path, files, par)  # 直接用训练好的就�?
            score_reads(path, files, contig_file, par, args.output)
            likelihood, read_probs = get_read_probs(path, files, args.output)
            df = pd.DataFrame(read_probs)
            df = df.T
            df[df < 1e-8] = 0
            hmm_file = os.path.join(os.path.abspath(os.path.dirname(args.output)), 'hmm_profile.tsv')
            df.to_csv(hmm_file, sep='\t', header=True)
        else:
            hmm_file = args.hmm_file
        X_hmm = gen_X_hmm(hmm_file, cov_file)

    bestK = gen_bestk(args.contig_file, X_t, args.clusters)
    logger.info("estimate clusters:\t" + str(bestK))

    logger.info("start calculate contig length")
    lengths = get_length(args.contig_file)
    length_weight = []
    for seq_id in namelist:
        length_weight.append(lengths[seq_id])

    seedURL = args.contig_file + ".seed"

    # 暂时不考虑整个contig_file 找不�?seed的极端情�?
    seed_list = []
    with open(seedURL) as f:
        for line in f:
            seed_list.append(line.rstrip('\n'))
    name_map = dict(zip(namelist, range(len(namelist))))
    # global seed_idx
    seed_idx = [name_map[seed_name] for seed_name in seed_list]


    # run weight kmeans
    logger.info("run kmeans with length weight")
    km = KMeans(n_clusters=bestK, init='k-means++', n_jobs=-1, n_init=30, random_state=7)
    km.fit(X_t, sample_weight=length_weight)
    idx = km.labels_
    kmeans_length_weight_output = os.path.dirname(args.output) + '/kmeans_length_weight_result.tsv'
    save_result(idx, kmeans_length_weight_output, namelist)
    kmeans_length_weight_output_dir = os.path.dirname(args.output) + '/kmeans_length_weight_result'
    os.mkdir(kmeans_length_weight_output_dir)
    gen_bins(contig_file, kmeans_length_weight_output, kmeans_length_weight_output_dir, "kmeans_weight_result")

    # run weight kmeans with composition information only
    logger.info("Run weight kmeans with composition information only.")
    km = KMeans(n_clusters=bestK, init='k-means++', n_jobs=-1, n_init=30, random_state=7)
    km.fit(X_com, sample_weight=length_weight)
    idx = km.labels_
    kmeans_length_weight_com_only_output = os.path.dirname(args.output) + '/kmeans_length_weight_com_only_result.tsv'
    save_result(idx, kmeans_length_weight_com_only_output, namelist)
    kmeans_length_weight_com_only_output_dir = os.path.dirname(args.output) + '/kmeans_length_weight_com_only_result'
    os.mkdir(kmeans_length_weight_com_only_output_dir)
    gen_bins(contig_file, kmeans_length_weight_com_only_output, kmeans_length_weight_com_only_output_dir, "com_result")

    # run weight kmeans with coverage information only
    if (len(X_cov[0]) >= 5):
        logger.info("Run weight kmeans with coverage information only.")
        km = KMeans(n_clusters=bestK, init='k-means++', n_jobs=-1, n_init=30, random_state=7)
        X_cov = np.log10(X_cov * int(100) + 1)  # 参考bisanity初始�?
        km.fit(X_cov, sample_weight=length_weight)
        idx = km.labels_
        kmeans_length_weight_cov_only_output = os.path.dirname(
            args.output) + '/kmeans_length_weight_cov_only_result.tsv'
        save_result(idx, kmeans_length_weight_cov_only_output, namelist)
        kmeans_length_weight_cov_only_output_dir = os.path.dirname(
            args.output) + '/kmeans_length_weight_cov_only_result'
        os.mkdir(kmeans_length_weight_cov_only_output_dir)
        gen_bins(contig_file, kmeans_length_weight_cov_only_output, kmeans_length_weight_cov_only_output_dir,
                 "cov_result")

    if args.pacbio_read_profiles:
        if (len(X_cov_pb[0]) >= 5):
            logger.info("Run weight kmeans partial seed with coverage pb information only.")
            km = KMeans(n_clusters=bestK, n_jobs=-1, n_init=30, random_state=7,
                            init=functools.partial(partial_seed_init, seed_idx=seed_idx))
            X_cov_pb_trans = np.log10(X_cov_pb * int(100) + 1)  # 参考bisanity初始�?
            km.fit(X_cov_pb, sample_weight=length_weight)
            idx = km.labels_
            kmeans_length_weight_cov_pb_only_partial_seed_output = os.path.dirname(
                args.output) + '/kmeans_length_weight_cov_pb_only_result_partial_seed.tsv'
            save_result(idx, kmeans_length_weight_cov_pb_only_partial_seed_output, namelist)
            kmeans_length_weight_cov_pb_only_partial_seed_output_dir = os.path.dirname(
                args.output) + '/kmeans_length_weight_cov_pb_only_result_partial_seed'
            os.mkdir(kmeans_length_weight_cov_pb_only_partial_seed_output_dir)
            gen_bins(contig_file, kmeans_length_weight_cov_pb_only_partial_seed_output,
                     kmeans_length_weight_cov_pb_only_partial_seed_output_dir, "cov_pb_result_partial_seed")

    if (len(X_cov_sr[0]) >= 5):
        logger.info("Run weight kmeans partial seed with coverage short read information only.")
        km = KMeans(n_clusters=bestK, n_jobs=-1, n_init=30, random_state=7,
                    init=functools.partial(partial_seed_init, seed_idx=seed_idx))
        X_cov_sr_trans = np.log10(X_cov_sr * int(100) + 1)  # 参考bisanity初始�?
        km.fit(X_cov_sr_trans, sample_weight=length_weight)
        idx = km.labels_
        kmeans_length_weight_cov_sr_only_partial_seed_output = os.path.dirname(
            args.output) + '/kmeans_length_weight_cov_sr_only_result_partial_seed.tsv'
        save_result(idx, kmeans_length_weight_cov_sr_only_partial_seed_output, namelist)
        kmeans_length_weight_cov_sr_only_partial_seed_output_dir = os.path.dirname(
            args.output) + '/kmeans_length_weight_cov_sr_only_result_partial_seed'
        os.mkdir(kmeans_length_weight_cov_sr_only_partial_seed_output_dir)
        gen_bins(contig_file, kmeans_length_weight_cov_sr_only_partial_seed_output,
                 kmeans_length_weight_cov_sr_only_partial_seed_output_dir, "cov_sr_result_partial_seed")

    # run kmeans with partial seed initial with length weight
    logger.info("Run kmeans sr with partial seed initial with length weight.")
    km = KMeans(n_clusters=bestK, n_jobs=-1, n_init=30, random_state=7,
                init=functools.partial(partial_seed_init, seed_idx=seed_idx))
    km.fit(np.hstack((X_cov_sr, X_com)), sample_weight=length_weight)
    idx = km.labels_
    kmeans_sr_partial_seed_length_weight_output = os.path.dirname(
        args.output) + '/kmeans_sr_partial_seed_length_weight_result.tsv'
    save_result(idx, kmeans_sr_partial_seed_length_weight_output, namelist)
    kmeans_sr_partial_seed_length_weight_output_dir = os.path.dirname(
        args.output) + '/kmeans_sr_partial_seed_length_weight_result'
    os.mkdir(kmeans_sr_partial_seed_length_weight_output_dir)
    gen_bins(contig_file, kmeans_sr_partial_seed_length_weight_output, kmeans_sr_partial_seed_length_weight_output_dir,
             "sr_partial_seed_length_weight_result")

    # run kmeans with partial seed initial with length weight
    if args.pacbio_read_profiles:
        logger.info("Run kmeans pb with partial seed initial with length weight.")
        km = KMeans(n_clusters=bestK, n_jobs=-1, n_init=30, random_state=7,
                    init=functools.partial(partial_seed_init, seed_idx=seed_idx))
        km.fit(np.hstack((X_cov_pb, X_com)), sample_weight=length_weight)
        idx = km.labels_
        kmeans_pb_partial_seed_length_weight_output = os.path.dirname(
            args.output) + '/kmeans_pb_partial_seed_length_weight_result.tsv'
        save_result(idx, kmeans_pb_partial_seed_length_weight_output, namelist)
        kmeans_pb_partial_seed_length_weight_output_dir = os.path.dirname(
            args.output) + '/kmeans_pb_partial_seed_length_weight_result'
        os.mkdir(kmeans_pb_partial_seed_length_weight_output_dir)
        gen_bins(contig_file, kmeans_pb_partial_seed_length_weight_output,
                 kmeans_pb_partial_seed_length_weight_output_dir, "pb_partial_seed_length_weight_result")

    # run weight kmeans with kmer_cov_hmm profile
    if args.use_hmm:
        logger.info("Run weight kmeans with kmer_cov_hmm information.")
        # nClass = sum((np.sum(X_hmm, axis=0) >= 100).astype(int)) #100 can be reset.
        km = KMeans(n_clusters=bestK, init='k-means++', n_jobs=-1, n_init=30, random_state=7)
        X_cov_kmer_hmm = np.hstack((X_t, X_hmm))
        km.fit(X_cov_kmer_hmm, sample_weight=length_weight)
        idx = km.labels_
        kmeans_length_weight_hmm_cov_kmer_output = os.path.dirname(
            args.output) + '/kmeans_length_weight_hmm_cov_kmer_result.tsv'
        save_result(idx, kmeans_length_weight_hmm_cov_kmer_output, namelist)
        kmeans_length_weight_hmm_cov_kmer_output_dir = os.path.dirname(
            args.output) + '/kmeans_length_weight_hmm_cov_kmer_result'
        os.mkdir(kmeans_length_weight_hmm_cov_kmer_output_dir)
        gen_bins(contig_file, kmeans_length_weight_hmm_cov_kmer_output, kmeans_length_weight_hmm_cov_kmer_output_dir,
                 "hmm_cov_kmer_weight_length_result")

    # run kmeans with partial seed initial with length weight
    logger.info("Run kmeans with partial seed initial with length weight.")
    km = KMeans(n_clusters=bestK, n_jobs=-1, n_init=30, random_state=7,
                init=functools.partial(partial_seed_init, seed_idx=seed_idx))
    km.fit(X_t, sample_weight=length_weight)
    idx = km.labels_
    kmeans_partial_seed_length_weight_output = os.path.dirname(
        args.output) + '/kmeans_partial_seed_length_weight_result.tsv'
    save_result(idx, kmeans_partial_seed_length_weight_output, namelist)
    kmeans_partial_seed_length_weight_output_dir = os.path.dirname(
        args.output) + '/kmeans_partial_seed_length_weight_result'
    os.mkdir(kmeans_partial_seed_length_weight_output_dir)
    gen_bins(contig_file, kmeans_partial_seed_length_weight_output, kmeans_partial_seed_length_weight_output_dir,
             "partial_seed_length_weight_result")

    # run kmeans with partial seed initial with length weight
    logger.info("Run kmeans with partial seed initial with length weight transform.")
    km = KMeans(n_clusters=bestK, n_jobs=-1, n_init=30, random_state=7,
                init=functools.partial(partial_seed_init, seed_idx=seed_idx))
    km.fit(np.log10(X_t * int(100) + 1), sample_weight=length_weight)
    idx = km.labels_
    kmeans_partial_seed_length_weight_output_trans = os.path.dirname(
        args.output) + '/kmeans_partial_seed_length_weight_transform_result.tsv'
    save_result(idx, kmeans_partial_seed_length_weight_output_trans, namelist)
    kmeans_partial_seed_length_weight_trans_output_dir = os.path.dirname(
        args.output) + '/kmeans_partial_seed_length_weight_transform_result'
    os.mkdir(kmeans_partial_seed_length_weight_trans_output_dir)
    gen_bins(contig_file, kmeans_partial_seed_length_weight_output_trans,
             kmeans_partial_seed_length_weight_trans_output_dir, "partial_seed_length_weight_trans_result")

    # run kmeans with partial seed initial
    logger.info("Run kmeans with partial seed initial.")
    km = KMeans(n_clusters=bestK, n_jobs=-1, n_init=30, random_state=7,
                init=functools.partial(partial_seed_init, seed_idx=seed_idx))
    km.fit(X_t)
    idx = km.labels_
    kmeans_partial_seed_initial_output = os.path.dirname(
        args.output) + '/kmeans_partial_seed_initial_ori_result.tsv'
    save_result(idx, kmeans_partial_seed_initial_output, namelist)
    kmeans_partial_seed_initial_output_dir = os.path.dirname(
        args.output) + '/kmeans_partial_seed_initial_ori_result'
    os.mkdir(kmeans_partial_seed_initial_output_dir)
    gen_bins(contig_file, kmeans_partial_seed_initial_output, kmeans_partial_seed_initial_output_dir,
             "partial_seed_initial_ori_result")

    checkm_out_dir = kmeans_partial_seed_initial_output_dir + '/checkm_out'
    os.mkdir(checkm_out_dir)
    output = args.output

    checkm_file = checkm_out_dir + "/checkm_analysis_init.txt"
    checkm_file_output = open((checkm_file), "w")
    threads = 45
    subprocess.call(["checkm", "lineage_wf", "-x", "bin", "-t", str(threads), kmeans_partial_seed_initial_output_dir,
                     checkm_out_dir], stdout=checkm_file_output)

    suffix_str = '.bin'
    checkm_analysis(checkm_file, suffix_str, checkm_out_dir)
    #
    goodbin_path = kmeans_partial_seed_initial_output_dir + '/good_bins'

    # 处理high_com_p_high_cont文件
    logger.info("Recluster the contigs from high_com_p_high_cont bins")
    high_com_p_high_cont_path = kmeans_partial_seed_initial_output_dir + "/High_completion_high_contamination"
    recluster_hh_bins(high_com_p_high_cont_path, mapObj, X_t, length_weight, namelist)

    # recluster other contigs
    logger.info("Recluster other contigs.")
    not_clustered_path = kmeans_partial_seed_initial_output_dir + "/others"
    recluster_other_contigs(not_clustered_path, X_t, namelist, mapObj, length_weight)

    convert(kmeans_partial_seed_initial_output_dir + '/good_bins', os.path.dirname(
        args.output) + '/kmeans_seed_partial_ori_with_postprocess.tsv')

    kmeans_seed_partial_with_postprocess_output_dir = os.path.dirname(
        args.output) + '/kmeans_seed_partial_with_postprocess_result'
    os.mkdir(kmeans_seed_partial_with_postprocess_output_dir)
    kmeans_seed_partial_with_postprocess_output = os.path.dirname(
        args.output) + '/kmeans_seed_partial_ori_with_postprocess.tsv'
    gen_bins(args.contig_file, kmeans_seed_partial_with_postprocess_output,
             kmeans_seed_partial_with_postprocess_output_dir, "seed_partial_with_postprocess")



    # run weight kmeans with kmer_cov_hmm profile _partial_seed
    if args.use_hmm:
        logger.info("Run weight kmeans partial seed with kmer_cov_hmm information.")
        # nClass = sum((np.sum(X_hmm, axis=0) >= 100).astype(int)) #100 can be reset.
        km = KMeans(n_clusters=bestK, n_init=30, random_state=7)
        X_cov_kmer_hmm = np.hstack((X_t, X_hmm))
        km.fit(X_cov_kmer_hmm, sample_weight=length_weight)
        idx = km.labels_
        kmeans_length_weight_hmm_cov_kmer_partial_seed_output = os.path.dirname(
            args.output) + '/kmeans_length_weight_hmm_cov_kmer_result_partial_seed.tsv'
        save_result(idx, kmeans_length_weight_hmm_cov_kmer_partial_seed_output, namelist)
        kmeans_length_weight_hmm_cov_kmer_partial_seed_output_dir = os.path.dirname(
            args.output) + '/kmeans_length_weight_hmm_cov_kmer_result_partial_seed'
        os.mkdir(kmeans_length_weight_hmm_cov_kmer_partial_seed_output_dir)
        gen_bins(contig_file, kmeans_length_weight_hmm_cov_kmer_partial_seed_output,
                 kmeans_length_weight_hmm_cov_kmer_partial_seed_output_dir,
                 "hmm_cov_kmer_weight_length_result_partial_seed")

        logger.info("Run weight kmeans partial seed with kmer_cov_hmm transform information.")
        # nClass = sum((np.sum(X_hmm, axis=0) >= 100).astype(int)) #100 can be reset.
        km = KMeans(n_clusters=bestK, n_jobs=-1, n_init=30, random_state=7,
                    init=functools.partial(partial_seed_init, seed_idx=seed_idx))
        X_cov_kmer_hmm = np.hstack((X_t, X_hmm))
        km.fit(np.log10(X_cov_kmer_hmm * int(100) + 1), sample_weight=length_weight)
        idx = km.labels_
        kmeans_length_weight_hmm_cov_kmer_partial_seed_trans_output = os.path.dirname(
            args.output) + '/kmeans_length_weight_hmm_cov_kmer_result_partial_seed_trans.tsv'
        save_result(idx, kmeans_length_weight_hmm_cov_kmer_partial_seed_trans_output, namelist)
        kmeans_length_weight_hmm_cov_kmer_partial_seed_trans_output_dir = os.path.dirname(
            args.output) + '/kmeans_length_weight_hmm_cov_kmer_result_partial_seed_trans'
        os.mkdir(kmeans_length_weight_hmm_cov_kmer_partial_seed_trans_output_dir)
        gen_bins(contig_file, kmeans_length_weight_hmm_cov_kmer_partial_seed_trans_output,
                 kmeans_length_weight_hmm_cov_kmer_partial_seed_trans_output_dir,
                 "hmm_cov_kmer_weight_length_result_partial_seed_trans")

    # run kmeans par
    seed_number_bacar_marker = gen_seed_number_bacar_marker(contig_file)
    if seed_number_bacar_marker > 0:
        logger.info("Run kmeans with partial bacar_marker seed initial with length weight transform.")
        seedURL = args.contig_file + ".bacar_marker.seed"
        seed_list = []
        with open(seedURL) as f:
            for line in f:
                seed_list.append(line.rstrip('\n'))
        name_map = dict(zip(namelist, range(len(namelist))))
        # global seed_idx
        seed_idx = [name_map[seed_name] for seed_name in seed_list]

        km = KMeans(n_clusters=seed_number_bacar_marker, n_jobs=-1, n_init=30, random_state=7,
                    init=functools.partial(partial_seed_init, seed_idx=seed_idx))
        km.fit(np.log10(X_t * int(100) + 1), sample_weight=length_weight)
        idx = km.labels_
        kmeans_partial_bacar_marker_seed_length_weight_output_trans = os.path.dirname(
            args.output) + '/kmeans_partial_bacar_marker_seed_length_weight_transform_result.tsv'
        save_result(idx, kmeans_partial_bacar_marker_seed_length_weight_output_trans, namelist)
        kmeans_partial_bacar_marker_seed_length_weight_trans_output_dir = os.path.dirname(
            args.output) + '/kmeans_partial_bacar_marker_seed_length_weight_transform_result'
        os.mkdir(kmeans_partial_bacar_marker_seed_length_weight_trans_output_dir)
        gen_bins(contig_file, kmeans_partial_bacar_marker_seed_length_weight_output_trans,
                 kmeans_partial_bacar_marker_seed_length_weight_trans_output_dir,
                 "partial_bacar_marker_seed_length_weight_trans_result")

    if seed_number_bacar_marker > 0:
        logger.info("Run kmeans with partial bacar_marker seed initial with length weight.")

        km = KMeans(n_clusters=seed_number_bacar_marker, n_jobs=-1, n_init=30, random_state=7,
                    init=functools.partial(partial_seed_init, seed_idx=seed_idx))
        km.fit(X_t, sample_weight=length_weight)
        idx = km.labels_
        kmeans_partial_bacar_marker_seed_length_weight_output = os.path.dirname(
            args.output) + '/kmeans_partial_bacar_marker_seed_length_weight_result.tsv'
        save_result(idx, kmeans_partial_bacar_marker_seed_length_weight_output, namelist)
        kmeans_partial_bacar_marker_seed_length_weight_output_dir = os.path.dirname(
            args.output) + '/kmeans_partial_bacar_marker_seed_length_weight_result'
        os.mkdir(kmeans_partial_bacar_marker_seed_length_weight_output_dir)
        gen_bins(contig_file, kmeans_partial_bacar_marker_seed_length_weight_output,
                 kmeans_partial_bacar_marker_seed_length_weight_output_dir,
                 "partial_bacar_marker_seed_length_weight_result")

    seedURL = args.contig_file + ".seed"
    # 暂时不考虑整个contig_file 找不�?seed的极端情�?
    seed_list = []
    with open(seedURL) as f:
        for line in f:
            seed_list.append(line.rstrip('\n'))
    name_map = dict(zip(namelist, range(len(namelist))))
    # global seed_idx
    seed_idx = [name_map[seed_name] for seed_name in seed_list]

    # das_tool 后处�?
    das_tool_dir = os.path.dirname(args.output) + "/das_tool_output"
    os.mkdir(das_tool_dir)
    if args.pacbio_read_profiles:
        if args.use_hmm:
            das_toolCmd = (
                        "DAS_Tool -i " + kmeans_length_weight_output + "," + kmeans_length_weight_com_only_output + ","
                        + kmeans_length_weight_cov_only_output + "," + kmeans_partial_seed_length_weight_output + ","
                        + kmeans_length_weight_hmm_cov_kmer_partial_seed_output + "," + kmeans_seed_partial_with_postprocess_output + "," + kmeans_partial_seed_length_weight_output_trans + "," + kmeans_pb_partial_seed_length_weight_output + "," + kmeans_sr_partial_seed_length_weight_output + "," + kmeans_length_weight_cov_sr_only_partial_seed_output + "," + kmeans_length_weight_cov_pb_only_partial_seed_output + "," + kmeans_partial_bacar_marker_seed_length_weight_output_trans + "," + kmeans_partial_bacar_marker_seed_length_weight_output + "," + kmeans_length_weight_hmm_cov_kmer_partial_seed_trans_output +
                        " -l kmeans_length_weight_result,kmeans_length_weight_com_only_result,kmeans_length_weight_cov_only_result,kmeans_partial_seed_length_weight_result,kmeans_length_weight_hmm_cov_kmer_result_partial_seed,kmeans_seed_partial_ori_with_postprocess,kmeans_partial_seed_length_weight_output_trans,kmeans_pb_partial_seed_length_weight_output,kmeans_sr_partial_seed_length_weight_output,kmeans_length_weight_cov_sr_only_partial_seed_output,kmeans_length_weight_cov_pb_only_partial_seed_output,kmeans_partial_bacar_marker_seed_length_weight_output_trans,kmeans_partial_bacar_marker_seed_length_weight_output,kmeans_length_weight_hmm_cov_kmer_partial_seed_trans_output"
                        + " -c " + args.contig_file + " -o " + das_tool_dir + "/das_tool_goodbins --threads 20 --score_threshold " + str(
                    args.binscore))
        else:
            das_toolCmd = (
                        "DAS_Tool -i " + kmeans_length_weight_output + "," + kmeans_length_weight_com_only_output + ","
                        + kmeans_length_weight_cov_only_output + "," + kmeans_partial_seed_length_weight_output + ","
                        + kmeans_seed_partial_with_postprocess_output + "," + kmeans_partial_seed_length_weight_output_trans + "," + kmeans_pb_partial_seed_length_weight_output + "," + kmeans_sr_partial_seed_length_weight_output + "," + kmeans_length_weight_cov_sr_only_partial_seed_output + "," + kmeans_length_weight_cov_pb_only_partial_seed_output + "," + kmeans_partial_bacar_marker_seed_length_weight_output_trans + "," + kmeans_partial_bacar_marker_seed_length_weight_output +
                        " -l kmeans_length_weight_result,kmeans_length_weight_com_only_result,kmeans_length_weight_cov_only_result,kmeans_partial_seed_length_weight_result,kmeans_seed_partial_ori_with_postprocess,kmeans_partial_seed_length_weight_output_trans,kmeans_pb_partial_seed_length_weight_output,kmeans_sr_partial_seed_length_weight_output,kmeans_length_weight_cov_sr_only_partial_seed_output,kmeans_length_weight_cov_pb_only_partial_seed_output,kmeans_partial_bacar_marker_seed_length_weight_output_trans,kmeans_partial_bacar_marker_seed_length_weight_output"
                        + " -c " + args.contig_file + " -o " + das_tool_dir + "/das_tool_goodbins --threads 20 --score_threshold " + str(
                    args.binscore))
    else:
        if args.use_hmm:
            das_toolCmd = (
                    "DAS_Tool -i " + kmeans_length_weight_output + "," + kmeans_length_weight_com_only_output + ","
                    + kmeans_length_weight_cov_only_output + "," + kmeans_partial_seed_length_weight_output + ","
                    + kmeans_length_weight_hmm_cov_kmer_partial_seed_output + "," + kmeans_seed_partial_with_postprocess_output + "," + kmeans_partial_seed_length_weight_output_trans + "," + kmeans_sr_partial_seed_length_weight_output + "," + kmeans_length_weight_cov_sr_only_partial_seed_output + "," + kmeans_partial_bacar_marker_seed_length_weight_output_trans + "," + kmeans_partial_bacar_marker_seed_length_weight_output + "," + kmeans_length_weight_hmm_cov_kmer_partial_seed_trans_output +
                    " -l kmeans_length_weight_result,kmeans_length_weight_com_only_result,kmeans_length_weight_cov_only_result,kmeans_partial_seed_length_weight_result,kmeans_length_weight_hmm_cov_kmer_result_partial_seed,kmeans_seed_partial_ori_with_postprocess,kmeans_partial_seed_length_weight_output_trans,kmeans_sr_partial_seed_length_weight_output,kmeans_length_weight_cov_sr_only_partial_seed_output,kmeans_partial_bacar_marker_seed_length_weight_output_trans,kmeans_partial_bacar_marker_seed_length_weight_output,kmeans_length_weight_hmm_cov_kmer_partial_seed_trans_output"
                    + " -c " + args.contig_file + " -o " + das_tool_dir + "/das_tool_goodbins --threads 20 --score_threshold " + str(
                args.binscore))
        else:
            das_toolCmd = (
                    "DAS_Tool -i " + kmeans_length_weight_output + "," + kmeans_length_weight_com_only_output + ","
                    + kmeans_length_weight_cov_only_output + "," + kmeans_partial_seed_length_weight_output + ","
                    + kmeans_seed_partial_with_postprocess_output + "," + kmeans_partial_seed_length_weight_output_trans + "," + kmeans_sr_partial_seed_length_weight_output + "," + kmeans_length_weight_cov_sr_only_partial_seed_output + "," + kmeans_partial_bacar_marker_seed_length_weight_output_trans + "," + kmeans_partial_bacar_marker_seed_length_weight_output +
                    " -l kmeans_length_weight_result,kmeans_length_weight_com_only_result,kmeans_length_weight_cov_only_result,kmeans_partial_seed_length_weight_result,kmeans_seed_partial_ori_with_postprocess,kmeans_partial_seed_length_weight_output_trans,kmeans_sr_partial_seed_length_weight_output,kmeans_length_weight_cov_sr_only_partial_seed_output,kmeans_partial_bacar_marker_seed_length_weight_output_trans,kmeans_partial_bacar_marker_seed_length_weight_output"
                    + " -c " + args.contig_file + " -o " + das_tool_dir + "/das_tool_goodbins --threads 20 --score_threshold " + str(
                args.binscore))

    logger.info("exec cmd: " + das_toolCmd)
    os.system(das_toolCmd)

    ###############

    # recluster remained contigs with kmeans and dbscan
    das_tool_output = das_tool_dir + '/das_tool_goodbins_DASTool_scaffolds2bin.txt'
    # recluster_contigs_after_das_tool(das_tool_output, namelist, mapObj, X_t, length_weight, l1_distance_score)
    clustered_contig_id = []
    with open(das_tool_output) as f:
        for line in f:
            clustered_contig_id.append(line.rstrip('\n').split('\t')[0])

    remained_contig_id = [namelist[i] for i in range(len(namelist)) if (namelist[i] not in clustered_contig_id)]

    outputdir = os.path.dirname(das_tool_output)
    prefix_str = 'remained_contigs_after_dastool'
    gen_remained_fasta_file(contig_file, remained_contig_id, outputdir, prefix_str)

    remained_contig_file = outputdir + "/" + prefix_str + '_0.bin'
    remained_contig_id_number = [mapObj[x] for x in remained_contig_id]
    X_t_remained = X_t[remained_contig_id_number]

    bin_number = gen_bestk(remained_contig_file, X_t_remained, 0)

    remained_contigs_weight = []
    for i in range(len(remained_contig_id_number)):
        remained_contigs_weight.append(length_weight[remained_contig_id_number[i]])

    seedURL = remained_contig_file + ".seed"
    global seed_idx
    if os.path.exists(seedURL):
        seed_list = []
        with open(seedURL) as f:
            for line in f:
                seed_list.append(line.rstrip('\n'))
        name_map = dict(zip(remained_contig_id, range(len(remained_contig_id))))
        seed_idx = [name_map[seed_name] for seed_name in seed_list]
        km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7,
                    init=functools.partial(partial_seed_init, seed_idx=seed_idx))
    else:
        km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7)

    km.fit(X_t_remained, sample_weight=remained_contigs_weight)
    idx = km.labels_
    save_result_refine(idx, remained_contig_file + ".reclustered.tsv",
                       namelist, remained_contig_id_number)


    output_dir = os.path.dirname(remained_contig_file) + '/remained_contig_file_result'
    os.mkdir(output_dir)
    gen_bins(contig_file, remained_contig_file + ".reclustered.tsv", output_dir, "remained_contig_file_result_result")

    checkm_out_dir = output_dir + '/checkm_out'
    os.mkdir(checkm_out_dir)
    output = args.output

    checkm_file = checkm_out_dir + "/checkm_analysis_init.txt"
    checkm_file_output = open((checkm_file), "w")
    threads = 20
    subprocess.call(["checkm", "lineage_wf", "-x", "bin", "-t", str(threads), output_dir,
                     checkm_out_dir], stdout=checkm_file_output)

    suffix_str = '.bin'
    checkm_analysis(checkm_file, suffix_str, checkm_out_dir)

    # deal with high_com_p_high_cont files
    logger.info("Recluster the contigs from high_com_p_high_cont bins")
    high_com_p_high_cont_path = output_dir + "/High_completion_high_contamination"
    recluster_hh_bins(high_com_p_high_cont_path, mapObj, X_t, length_weight, namelist)

    # recluster other contigs
    logger.info("Recluster other contigs.")
    not_clustered_path = output_dir + "/others"
    recluster_other_contigs(not_clustered_path, X_t, namelist, mapObj, length_weight)

    convert(output_dir + '/good_bins', os.path.dirname(
        args.output) + '/remianed_contigs_postprocess.tsv')

    f = open(args.output, 'w')
    f.write("@Version:0.9.0\n")
    f.write("@SampleID:{}\n".format(os.path.basename(contig_file)))
    f.write("@@SEQUENCEID\tBINID\n")

    file_1 = open(das_tool_output, 'r')
    file_2 = open(os.path.dirname(
        args.output) + '/remianed_contigs_postprocess.tsv', 'r')

    lines_1 = file_1.readlines()
    file_1.close()
    lines_2 = file_2.readlines()
    file_2.close()

    f.writelines(lines_1)
    f.writelines(lines_2)
    f.close()
