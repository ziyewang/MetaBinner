#!/usr/bin/env python

#modified from https://github.com/dparks1134/UniteM/blob/master/unitem/ensemble.py
import copy
import os
import random
import itertools
from collections import defaultdict
import sys

from metabinner_util import get_bin_dirs

from unitem_common import (read_bins,
                           parse_bin_stats,
                           calculateN50L50M50)

from biolib.common import (make_sure_path_exists,
                           check_dir_exists,
                           check_file_exists,
                           query_yes_no,
                           remove_extension)

from unitem_markers import Markers

#
em_mode = 'greedy'  # sys.argv[1] #'combo'
bac_mg_table = sys.argv[1]
ar_mg_table = sys.argv[2]
bin_dirs_file = sys.argv[3]
output_dir = sys.argv[4]
####
# parameters

CONT_WEIGHT = 3  # help="weight given to contamination for assessing genome quality", type=float, default=2)
SEL_MIN_QUALITY = 50.0  # add_argument('-q', '--sel_min_quality', help="minimum quality of bin to consider during bin selection process", type=float, default=50)
SEL_MIN_COMP = 50.0  # add_argument('-x', '--sel_min_comp', help="minimum completeness of bin to consider during bin selection process", type=float, default=50)
SEL_MAX_CONT = 15.0  # add_argument('-y', '--sel_max_cont', help="maximum contamination of bin to consider during bin selection process", type=float, default=10)
REMOVE_PERC = 50.0  # add_argument('-r', '--remove_perc', help="minimum percentage of bins required to remove contigs from highest-quality bin", type=float, default=50.0)
ADD_PERC = 50.0  # add_argument('-a', '--add_perc', help="minimum percentage of matched bins required to add contigs to highest-quality bin", type=float, default=50.0)
ADD_MATCHES = 2  # add_argument('-m', '--add_matches', help="minimum number of matched bins required to 'add' contigs", type=int, default=2)

REPORT_MIN_QUALITY = 10.0  # add_argument('--report_min_quality', help="minimum quality of bin to report", type=float, default=10)

sel_q = 'sel_q_' + str(SEL_MIN_QUALITY)

cont_weight = '_cont_weight_' + str(CONT_WEIGHT)
reprot_minq = '_reprot_minq_' + str(REPORT_MIN_QUALITY)

add_matches = '_add_matches_' + str(ADD_MATCHES)

overlap_strategy = 'unitem_overlap_strategy'  # default

if overlap_strategy == 'unitem_overlap_strategy':
    MIN_PERC_COMMON = 0.5
elif overlap_strategy == 'metawrap_overlap_strategy':
    MIN_PERC_COMMON = 0.8

quality_filter = False  # default True
if quality_filter:
    quality_filter_flag = 'True'
else:
    quality_filter_flag = 'False_update'

# sort_quality_strategy = 'sort_matches_max_q'
# defalut='sort_matches_sum_q'
sort_quality_strategy = 'sort_matches_maxaddsum_q'

BIN_PREFIX = 'bin'


####################################################################################################
################################# different component methods          ############################
####################################################################################################

# megahit short read solidbin cov
# bin_dirs_file = '/home/wangzy/data/cami2019/cami2b/output/use_solidbin_scriptes_get_kmer_cov/rhimgCAMI2_short_read_pooled_megahit_map/output/metabinner_multiseed_update/component_result.X_t_logtrans_bin_dirs.tsv'
# output_dir='/home/wangzy/data/cami2019/cami2b/output/use_solidbin_scriptes_get_kmer_cov/rhimgCAMI2_short_read_pooled_megahit_map/output/new_ensemble_strategy/em_res_X_t_logtrans'

class MetabinnerEnsemble():
    """
      Bin selection:
    consensus -> Consensus clustering across multiple binning methods #
    greedy    -> Greedy bin selection across multiple binning methods （DAS_tool）
    unanimous -> Unanimous bin filtering across multiple binning methods (Binning_refiner )
    """

    def __init__(self):
        """

        """

    def _matched_bin_sets(self, bins, contig_lens, bin_quality, min_quality, min_comp, max_cont, no_bin_matching,
                          overlap_strategy='unitem_overlap_strategy', sort_quality_strategy='sort_matches_sum_q',
                          quality_filter=True):
        """Determine all sets of matched bins.

        Parameters
        ----------
        bins : d[binning method][bin ID]
          Contigs for bins across all binning methods.
        contig_lens : d[cid] -> length of contig
          Length of contigs.
        bin_quality : list of tuples
          Bin quality information. Must be sorted by quality!
        min_quality : float
          Minimum quality of bin to consider during bin matching.
        no_bin_matching : boolean
          Flag indicating bin matching should be skipped and all bins treated independently.

        Return
        ------
          Matched bin sets.
        """

        matched_sets = []
        processed_bins = defaultdict(set)
        for cur_bm, cur_bid, _domain, comp, cont, quality, N50, gs in bin_quality:
            if cur_bid in processed_bins[cur_bm]:
                continue  # bin has already been considered

            if not bins[cur_bm][cur_bid]:
                continue  # null bin
            matched_bins = []
            matched_bins.append((cur_bm, cur_bid, quality, N50, gs))
            if quality_filter:
                quality_flag = (quality >= min_quality and comp >= min_comp and cont <= max_cont)
            else:
                quality_flag = True
            if not no_bin_matching and quality_flag:
                binned_bm = []
                binned_bm.append(cur_bm)
                for test_bm, test_bid, _domain, comp, cont, quality, N50, gs in bin_quality:
                    if test_bm in binned_bm:
                        # can't group with a bin from the same binning method
                        continue

                    if not bins[test_bm][test_bid]:
                        continue  # null bin

                    if test_bid in processed_bins[test_bm]:
                        continue  # bin has already been considered

                    if quality_filter:
                        if quality < min_quality or comp < min_comp or cont > max_cont:
                            continue

                    bp_in_common, total_bp1, total_bp2 = self._bases_in_common(bins[cur_bm][cur_bid],
                                                                               bins[test_bm][test_bid],
                                                                               contig_lens)

                    if overlap_strategy == 'unitem_overlap_strategy':
                        per_bp_in_common = bp_in_common / max(total_bp1, total_bp2)
                        if per_bp_in_common > MIN_PERC_COMMON:
                            matched_bins.append((test_bm, test_bid, quality, N50, gs))
                            binned_bm.append(test_bm)
                    elif overlap_strategy == 'metawrap_overlap_strategy':
                        per_bp_in_common = bp_in_common / min(total_bp1, total_bp2)
                        if per_bp_in_common > MIN_PERC_COMMON:
                            matched_bins.append((test_bm, test_bid, quality, N50, gs))
                            binned_bm.append(test_bm)

            # removed matched bins
            for bm, bid, _q, _n50, _gs in matched_bins:
                processed_bins[bm].add(bid)

            matched_sets.append(tuple(matched_bins))

        # sort by total quality of bins in matched sets
        if sort_quality_strategy == 'sort_matches_max_q':
            matched_sets.sort(key=lambda ms: (max([x[2] for x in ms]),
                                              max([x[3] for x in ms]),
                                              max([x[4] for x in ms]),
                                              random.random()),
                              reverse=True)
        elif sort_quality_strategy == 'sort_matches_sum_q':
            matched_sets.sort(key=lambda ms: (sum([x[2] for x in ms]),
                                              sum([x[3] for x in ms]),
                                              sum([x[4] for x in ms]),
                                              random.random()),
                              reverse=True)
        elif sort_quality_strategy == 'sort_matches_maxaddsum_q':
            matched_sets.sort(key=lambda ms: (max([x[2] for x in ms]),
                                              sum([x[2] for x in ms]),
                                              max([x[3] for x in ms]),
                                              max([x[4] for x in ms]),
                                              random.random()),
                              reverse=True)
        return matched_sets

    def _bases_in_common(self, contig_ids1, contig_ids2, contig_lens):
        """Calculate number of base pairs in common between two bins.

        Parameters
        ----------
        contig_ids1 : iterable
          Contigs in bin.
        contig_ids2 : iterable
          Contigs in bin.
        contig_lens : d[seq ID] -> length
          Length of contigs.

        Returns
        -------
          Base pairs in common, total bases in bin 1, total bases in bin 2
        """

        bp_in_common = sum([contig_lens[seq_id]
                            for seq_id in contig_ids2
                            if seq_id in contig_ids1])

        total_bp1 = sum([contig_lens[seq_id] for seq_id in contig_ids1])
        total_bp2 = sum([contig_lens[seq_id] for seq_id in contig_ids2])

        return bp_in_common, total_bp1, total_bp2

    # from unitem Ensemble class
    # update for metabinner
    def _bin_quality(self, markers, bins, contigs, gene_tables, quality_weight):
        """Determine estimated completeness, contamination, and quality of bins.

        Parameters
        ----------
        bins : d[binning method][bin ID] -> set(cid1, cid2, ... cidN)
          Contigs for bins across all binning methods.
        contigs : d[cid] -> seq
          Contigs across all bins.
        gene_tables : d[binning method] -> (bac gene table, ar gene table)
          Bacterial and archaeal marker gene tables for bins in each binning method.
        quality_weight : float
          Weight given to contamination when assessing genome quality.

        Return
        ------
          List with bin metadata sorted by quality, then N50, the genome size.
        """

        # markers = Markers()

        q = []

        for method_id in methods_sorted:
            for bin_id in bins[method_id]:
                domain, comp, cont = markers.bin_quality(bins[method_id][bin_id])
                bin_seqs = [contigs[cid] for cid in bins[method_id][bin_id]]
                n50, l50, m50 = calculateN50L50M50(bin_seqs)
                genome_size = sum([len(s) for s in bin_seqs])
                q.append((method_id, bin_id, domain, comp, cont, comp - quality_weight * cont, n50, genome_size))

        # sort bins by quality follwed by N50 followed by genome size, and
        # break remaining ties randomly
        q.sort(key=lambda x: (x[5], x[6], x[7], random.random()),
               reverse=True)

        return q

    def _reconstruct_match_sets(self, markers, matched_set,
                                bins,
                                contigs,
                                bin_quality,
                                remove_perc,
                                add_perc,
                                add_matches,
                                sel_min_quality,
                                sel_min_comp,
                                sel_max_cont,
                                em_mode='greedy', report_min_quality=REPORT_MIN_QUALITY):
        """
        changed from unitem ensemble.py def _resolve_matched_set()
        Select contigs to cluster from matched bin set.

        Parameters
        ----------
        matched_set : iterable with bin metadata (binning method, bin ID, quality)
        Matched set of bins. Highest quality bin must be first.
        bins : d[binning method][bin ID]
        Contigs for bins across all binning methods.
        contigs : d[cid] -> seq
        Contigs across all bins.
        bin_quality : list of tuples
        Bin quality information. Must be sorted by quality!!!!!
        remove_perc : float
        Minimum percentage of bins from other binning methods require to remove contig in highest quality bin.
        add_perc : float
        Minimum percentage of matched bins required to add contig to highest quality bin.
        add_matches : float
        Minimum number of matched bins required to 'add' contigs.
        sel_min_quality : float
        Minimum quality of bins to consider for filtering contigs.
        """
        # identify contig count for bins in matched in set

        primary_bm, primary_bid, primary_q, primary_n50, primary_gs = matched_set[0]
        primary_contigs = bins[primary_bm][primary_bid]

        if primary_q == 100:
            new_bin = {}
            removed_contigs = set()
            added_contigs = set()
            for cid in primary_contigs:
                new_bin[cid] = contigs[cid]

        else:
            matched_contigs = defaultdict(int)
            matched_bin_ids = set()
            matched_methods = set()
            for bm, bid, _q, _n50, _gs in matched_set:
                matched_bin_ids.add(bm + bid)
                matched_methods.add(bm)
                for cid in bins[bm][bid]:
                    matched_contigs[cid] += 1

            new_bin = {}
            removed_contigs = set()
            added_contigs = set()

            if em_mode == 'greedy':
                for cid in primary_contigs:
                    new_bin[cid] = contigs[cid]

        return new_bin, removed_contigs, added_contigs

    def _update_bins(self, bins, cids_to_remove):
        """Remove specified contigs from all bins.

            Parameters
            ----------
            bins : d[binning method][bin ID]
              Contigs for bins across all binning methods.
            cids_to_remove : iterable
              Contigs to remove from bins.
        """

        cids_to_remove = set(cids_to_remove)

        # remove scaffolds in highest quality bin from all other bins
        for method_id in bins:
            for bin_id in bins[method_id]:
                bins[method_id][bin_id] -= cids_to_remove


# for each bin
# update for checkm marker
def get_binstats(bm, bin_id, bin_contig_names, length_dict):
    bin_seqs = [contigs[cid] for cid in bin_contig_names]
    n50, l50, m50 = calculateN50L50M50(bin_seqs)
    genome_size = sum([len(s) for s in bin_seqs])

    domain, comp, cont = markers.bin_quality(bin_contig_names)
    cc_score = comp - cont * CONT_WEIGHT

    return (bm, bin_id, domain, comp, cont, cc_score, n50, genome_size)


# update for checkm marker
def get_init_bin_quality(orig_bins, methods_sorted, contig_lens):
    # calculate quality of each bin of the orig_bins.
    init_bin_quality = defaultdict(lambda: {})
    for method_id in methods_sorted:
        for bin_id in orig_bins[method_id]:
            init_bin_quality[method_id][bin_id] = get_binstats(method_id, bin_id, orig_bins[method_id][bin_id],
                                                               contig_lens)

    return init_bin_quality


def write_init_bin_quality_per_bin(init_bin_quality, methods_sorted, output_dir):
    """
     #unitem：Record initial state of each contig；#converted to save per_bin's quality by me
    :param init_bin_quality:
    :param methods_sorted:
    :return:
    """
    fout = open(os.path.join(output_dir, 'bin_info_initial.tsv'), 'w')
    fout.write('Binning_method\tBin ID\tdomain\tcomp\tcont\tcc_score\tN50\tgenome_size\n')
    for method_id in methods_sorted:
        for bin_id in init_bin_quality[method_id]:
            fout.write('%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % init_bin_quality[method_id][bin_id])

    fout.close()


markers = Markers()

bin_dirs = get_bin_dirs(bin_dirs_file)

bins, contigs, contigs_in_bins = read_bins(bin_dirs)


methods_sorted = sorted(bins.keys())
contig_lens = {cid: len(contigs[cid]) for cid in contigs}
orig_bins = copy.deepcopy(bins)

gene_tables = markers.marker_gene_tables(bac_mg_table, ar_mg_table)

# create output directories


bin_dir = os.path.join(output_dir, em_mode + cont_weight + '_mincomp_' + str(SEL_MIN_COMP) + '_maxcont_' + str(
    SEL_MAX_CONT) + '_bins')
print("bin_dir:\t" + bin_dir)

make_sure_path_exists(bin_dir)

# get initial quality of bins
print('Get initial quality of bins.')
init_bin_quality = get_init_bin_quality(orig_bins, methods_sorted, contig_lens)

write_init_bin_quality_per_bin(init_bin_quality, methods_sorted, output_dir)

# ensemble of binning results
ensemble = MetabinnerEnsemble()
bin_num = 0
total_comp = 0
total_cont = 0
total_quality = 0
sel_gene_tables = {}
sel_bins = {}
selected_rows = {}
out_bin_quality = {}
######################################################################
########################## greedy  ##############################
######################################################################
below_min1_count = 0
no_bin_matching = True  # greedy

if em_mode != 'greedy':
    no_bin_matching = False  # not greedy

out_file = bin_dir + '_res.tsv'
with open(out_file, 'w') as handle_file:
    # step 1: greedy
    while True:
        # determine highest quality match bin set
        bins_quality = ensemble._bin_quality(markers, bins, contigs, gene_tables, quality_weight=CONT_WEIGHT)
        matched_sets = ensemble._matched_bin_sets(bins,
                                                  contig_lens,
                                                  bins_quality,
                                                  SEL_MIN_QUALITY,
                                                  SEL_MIN_COMP,
                                                  SEL_MAX_CONT,
                                                  no_bin_matching, overlap_strategy=overlap_strategy,
                                                  sort_quality_strategy=sort_quality_strategy,
                                                  quality_filter=quality_filter)

        if len(matched_sets) == 0:
            break  # no bins to be resolve

        new_bin, removed_contigs, added_contigs = ensemble._reconstruct_match_sets(markers, matched_sets[0],
                                                                                   bins,
                                                                                   contigs,
                                                                                   bins_quality,
                                                                                   REMOVE_PERC,
                                                                                   ADD_PERC,
                                                                                   ADD_MATCHES,
                                                                                   SEL_MIN_QUALITY,
                                                                                   SEL_MIN_COMP,
                                                                                   SEL_MAX_CONT,
                                                                                   em_mode=em_mode)

        domain, comp, cont = markers.bin_quality(new_bin.keys())
        cc_score = comp - cont * CONT_WEIGHT

        if cc_score < REPORT_MIN_QUALITY:
            break
        else:
            if (comp >= SEL_MIN_COMP) and (cont <= SEL_MAX_CONT):
                total_comp += comp
                total_cont += cont
                total_quality += cc_score

                # report selection
                bin_num += 1
                out_bin_in = '%s_%d' % (BIN_PREFIX, bin_num)
                out_bin_quality[out_bin_in] = (comp, cont)
                primary_bm, primary_bid, _q, _n50, _gs = matched_sets[0][0]
                print("Selected %s from %s with quality = %.1f (comp. = %.1f%%, cont. = %.1f%%)." % (
                    primary_bid,
                    primary_bm,
                    cc_score,
                    comp,
                    cont))

                if em_mode != 'greedy' and em_mode != 'unanimous':
                    # performing consensus binning
                    print("-> Identified %d matched bins, removed %d contigs, added %d contigs." % (
                        len(matched_sets[0]),
                        len(removed_contigs),
                        len(added_contigs)))
                elif em_mode != 'greedy':
                    # performing unanimous binning
                    print("-> Identified %d matched bins and removed %d contigs." % (
                        len(matched_sets[0]),
                        len(removed_contigs)))

                # write out contig info
                matched_bins = {}
                for m in matched_sets[0]:
                    matched_bins[m[0]] = m[1]

                new_bin_file = os.path.join(bin_dir, out_bin_in + '.fna')
                fout_bin = open(new_bin_file, 'w')
                for cid, seq in new_bin.items():
                    fout_bin.write('>%s\n' % cid)
                    fout_bin.write(seq + '\n')

                fout_bin.close()

                for key in new_bin.keys():
                    handle_file.write(key + '\t' + str(bin_num) + '\n')

            ensemble._update_bins(bins, new_bin.keys())
