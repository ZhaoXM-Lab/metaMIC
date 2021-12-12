#!/usr/bin/env python

import pandas as pd
import numpy as np
import multiprocessing
import argparse
import operator
import os
import random
import sys
import time
import random
import subprocess
import pysam
import collections
import warnings
import math
import re
from Bio import SeqIO

base_path = os.path.split(__file__)[0]

def fragment_distribution(samfile):
    all_reads = samfile.fetch()
    size_freq = collections.defaultdict(int)
    for read in all_reads:
        if read.rnext == read.tid and read.is_paired:
            size = abs(read.isize)
            size_freq[size] += 1
    return size_freq


def FragMAD(freq):
    """
    calculate median and median absolute deviation fragment size distribution
    """
    all_size = []
    for key, value in freq.items():
        all_size.extend([key] * int(value))
    median_size = np.median(all_size)
    residuals = abs(np.array(all_size) - median_size)
    mad_size = 1.4826 * np.median(residuals)
    return median_size, mad_size


def split_sam(args):
    split_command = ' '.join(['sh',
                              os.path.join(base_path, "split_sam.sh"),
                              args.assemblies,
                              args.bamfile,
                              args.output,
                              args.samtools])
    os.system(split_command)


def seq_parse(args):
    input = SeqIO.parse(args.assemblies, "fasta")
    contig_seqs = {}
    for record in input:
        if len(record.seq) >= args.min_length:
            contig_seqs[record.id] = str(record.seq)
    return contig_seqs


def kmer_parse(seq, pool):
    seq_kmer = {"position": [], "KAD": []}
    for i in range(len(seq)):
        if seq[i:(i + 25)] in pool:
            seq_kmer["KAD"].append(pool[seq[i:(i + 25)]])
            seq_kmer["position"].append(i + 1)
        if (i + 25) >= len(seq):
            break
    return seq_kmer


def KAD_window_cal(seq_kmer):
    KAD_window_dict = {"start_pos": [],
                       "mean_KAD": [],
                       "abnormal_KAD_ratio": [],
                       "dev_KAD": []}

    for i in range(300, len(seq_kmer['position']), 100):
        KAD_window_dict["start_pos"].append(i)
        mean_KAD = np.mean(np.abs(seq_kmer['KAD'][i:i + 100]))
        KAD_window_dict["mean_KAD"].append(mean_KAD)
        KAD_window_dict["abnormal_KAD_ratio"].append(
            np.sum(np.abs(seq_kmer['KAD'][i:i + 100]) > 0.5) / 100)
        KAD_window_dict["dev_KAD"].append(
            np.sqrt(np.var(np.abs(seq_kmer['KAD'][i:i + 100]))))
    return KAD_window_dict


def KAD_feature(args):
    seq_data = seq_parse(args)
    KAD_dict = {"contig": [],
                'start_pos': [],
                'mean_KAD': [],
                'abnormal_KAD_ratio': [],
                'dev_KAD': []}
    for contig, seq in seq_data.items():
        if len(seq) < args.min_length:
            continue
        if os.path.exists(os.path.join(args.output, "temp/KAD/KAD_data/",
                                                    "{}.KAD".format(str(contig)))):
            try:
                KAD_data = pd.read_csv(os.path.join(args.output, "temp/KAD/KAD_data/",
                                                    "{}.KAD".format(str(contig))), index_col=0, sep="\t")
                KAD_data = KAD_data.drop_duplicates(['k-mer'])
            except BaseException:
                continue
            KAD_data.index = KAD_data['k-mer']
            KAD_pool = KAD_data.loc[:, 'KAD'].to_dict()
            seq_kmer = kmer_parse(seq, KAD_pool)
            KAD_window = KAD_window_cal(seq_kmer)
            KAD_dict["contig"].extend([contig] * len(KAD_window['start_pos']))
            KAD_dict["start_pos"].extend(KAD_window['start_pos'])
            KAD_dict["mean_KAD"].extend(KAD_window["mean_KAD"])
            KAD_dict["abnormal_KAD_ratio"].extend(
                KAD_window["abnormal_KAD_ratio"])
            KAD_dict["dev_KAD"].extend(KAD_window["dev_KAD"])
    return KAD_dict


def KAD(args, contig, file):
    if os.path.exists(os.path.join(args.output, "temp/KAD/KAD_data/",
                                   str(contig), ".KAD")):
        return 0
    contig_file = os.path.join(args.output, "temp/split/contigs/", "{}.fa".format(file))

    read_file = os.path.join(args.output,
                             "temp/split/reads/{}.read.fa".format(str(contig)))
    # kmer count
    outputdir = os.path.join(args.output, "temp/KAD/temp")
    contig_command1 = ' '.join([args.jellyfish,
                                "count -m 25 -o",
                                os.path.join(outputdir, '{}.jf'.format(str(contig))),
                                "-s 100M -t 8",
                                contig_file])
    contig_command2 = ' '.join([args.jellyfish,
                                "dump -c -t -o",
                                os.path.join(outputdir, '{}_count.txt'.format(str(contig))),
                                os.path.join(outputdir, '{}.jf'.format(str(contig)))])
    os.system(contig_command1)
    os.system(contig_command2)
    read_command1 = ' '.join([args.jellyfish,
                              "count -m 25 -o",
                              os.path.join(outputdir, '{}.read.jf'.format(str(contig))),
                              "-s 100M -t 8",
                              read_file])
    read_command2 = ' '.join([args.jellyfish,
                              "dump -c -t -o",
                              os.path.join(outputdir, '{}_count.read.txt'.format(str(contig))),
                              os.path.join(outputdir, '{}.read.jf'.format(str(contig)))])
    os.system(read_command1)
    os.system(read_command2)
    assembly_kmer = pd.read_csv(os.path.join(args.output, "temp/KAD/temp/",
                                             "{}_count.txt".format(str(contig))), sep="\t", header=None)
    assembly_kmer.index = assembly_kmer[0]

    try:
        read_kmer = pd.read_csv(os.path.join(args.output, "temp/KAD/temp/",
                                             "{}_count.read.txt".format(str(contig))),
                                             sep="\t", header=None)
        read_kmer.index = read_kmer[0]
    except BaseException:
        # zero reads mapped to contig
        return 0
    shared_kmer = set(assembly_kmer.loc[assembly_kmer[1] == 1, 0]).intersection(read_kmer.index)

    if len(shared_kmer) == 0:
        kmer_depth = pd.value_counts(read_kmer.loc[read_kmer[1] > 5, 1]).index[0]
    else:
        kmer_depth = pd.value_counts(read_kmer.loc[shared_kmer, ][1]).index[0]

    assembly_kmer.columns = ['k-mer', 'assembly_count']
    read_kmer.columns = ['k-mer', 'read_count']
    assembly_kmer.index = range(assembly_kmer.shape[0])
    read_kmer.index = range(read_kmer.shape[0])
    kmer_result = pd.merge(assembly_kmer, read_kmer, how='outer')
    kmer_result = kmer_result.fillna(0)
    kmer_result['KAD'] = np.log2((kmer_result['read_count'] + kmer_depth)
                                 / (kmer_depth * (kmer_result['assembly_count'] + 1)))
    kmer_result.loc[(kmer_result['read_count'] == 1) *
                    (kmer_result['assembly_count'] == 0), 'KAD'] = np.nan
    kmer_result = kmer_result.loc[kmer_result['KAD'] == kmer_result['KAD'], ]
    kmer_result.loc[:, ['k-mer', 'KAD']].to_csv(
        os.path.join(args.output, "temp/KAD/KAD_data/", "{}.KAD".format(str(contig))), sep="\t")


def fragment_coverage_cal(reads, mu, dev, length):
    """
    calculate fragment coverage per contig
    """
    frag_coverage = np.array([0] * length)
    for read in reads:
        if read.rnext == read.tid and read.is_proper_pair:
            size = abs(read.isize)
            if (mu - 3 * dev <= size <= mu + 3 * dev):
                if read.next_reference_start < read.reference_start:
                    start = min(read.next_reference_start,
                                read.reference_start,
                                read.reference_end)
                    end = start + size
                    frag_coverage[start:end] += 1
    return frag_coverage


def window_read_cal(reads, mu, dev):
    read_dict = {"start_pos": [], "read_count": [], "proper_read_count": [], "inversion_read_count": [], "clipped_read_count": [],
                 "supplementary_read_count": [], "discordant_size_count": [], "discordant_loc_count": []}
    read_temp = {"num_read": 0, "num_proper": 0, "num_inversion": 0, "num_clipped": 0, "num_supplementary": 0, "num_discordant_size": 0,
                 "num_discordant_loc": 0}
    pos = 0
    for read in reads:
        new_pos = math.floor((read.reference_start - 300) / 100) * 100 + 300
        if read.reference_start < 300:
            continue
        if pos == 0:
            pos = new_pos
        elif new_pos != pos:
            read_dict["start_pos"].append(pos)
            read_dict["read_count"].append(read_temp["num_read"])
            read_dict["proper_read_count"].append(read_temp["num_proper"])
            read_dict["inversion_read_count"].append(
                read_temp["num_inversion"])
            read_dict["clipped_read_count"].append(read_temp["num_clipped"])
            read_dict["supplementary_read_count"].append(
                read_temp["num_supplementary"])
            read_dict["discordant_size_count"].append(
                read_temp["num_discordant_size"])
            read_dict["discordant_loc_count"].append(
                read_temp["num_discordant_loc"])
            read_temp = {"num_read": 0,
                         "num_proper": 0,
                         "num_inversion": 0,
                         "num_clipped": 0,
                         "num_supplementary": 0,
                         "num_discordant_size": 0,
                         "num_discordant_loc": 0}
            pos = new_pos
        read_temp["num_read"] += 1
        if read.is_paired:
            if read.rnext == read.tid:
                if read.is_proper_pair:
                    read_temp["num_proper"] += 1
                if (read.is_reverse + read.mate_is_reverse) != 1:
                    read_temp["num_inversion"] += 1
                if not mu - 3 * dev <= abs(read.isize) <= mu + 3 * dev:
                    read_temp["num_discordant_size"] += 1
            else:
                read_temp["num_discordant_loc"] += 1
        if read.get_cigar_stats()[0][4] > 20:
            read_temp["num_clipped"] += 1
        if (read.is_supplementary and read.get_cigar_stats()[0][5] > 20):
            read_temp["num_supplementary"] += 1
    return read_dict


def window_frag_cal(coverage):
    """
    Using sliding window approach to smooth out features
    """
    coverage = np.array(coverage)
    cov = {"pos": [], "coverage": [], "deviation": []}
    for i in range(300, len(coverage), 100):
        start = i
        end = i + 100
        cov["coverage"].append(np.mean(coverage[start:end]))
        cov["deviation"].append(
            np.sqrt(np.var(coverage[start:end])) / np.mean(coverage[start:end]))
        cov["pos"].append(start)
        if len(coverage) - end <= 300:
            break
    return cov


def contig_pool(samfile):
    contig_len = {}
    for (ref, lens) in zip(samfile.references, samfile.lengths):
        contig_len[ref] = lens
    return contig_len


def pileup_window_cal(pileup_dict):
    window_dict = {"contig": [], "start_pos": [], "correct_portion": [], "ambiguous_portion": [], "disagree_portion": [],
                   "deletion_portion": [], "insert_portion": [], "coverage": [], "deviation": []}
    for i in range(300, len(pileup_dict['correct']), 100):
        start = i
        end = i + 100
        total = np.sum(pileup_dict['depth'][start:end])
        window_dict["contig"].append(pileup_dict["contig"][0])
        window_dict["start_pos"].append(start)
        window_dict["correct_portion"].append(
            np.sum(pileup_dict['correct'][start:end]) / total)
        window_dict["ambiguous_portion"].append(
            np.sum(pileup_dict["ambiguous"][start:end]) / total)
        window_dict["insert_portion"].append(
            np.sum(pileup_dict['insert'][start:end]) / total)
        window_dict["deletion_portion"].append(
            np.sum(pileup_dict['deletion'][start:end]) / total)
        window_dict["disagree_portion"].append(
            np.sum(pileup_dict['disagree'][start:end]) / total)
        window_dict["coverage"].append(
            np.mean(pileup_dict["depth"][start:end]))
        window_dict["deviation"].append(np.sqrt(np.var(
            pileup_dict["depth"][start:end])) / np.mean(pileup_dict["depth"][start:end]))
        if len(pileup_dict['correct']) - (i + 100) <= 300:
            break
    return window_dict


def read_breakpoint_per_contig(samfile, ref, lens):
    reads = samfile.fetch(contig=ref)
    break_count = {"breakcount": np.array([0] * lens),
                   "readcount": np.array( [0] * lens)}
    for read in reads:
        ref_end = read.reference_end
        ref_start = read.reference_start
        read_start = read.query_alignment_start
        read_end = read.query_alignment_end
        break_count["readcount"][ref_start:ref_end] += 1

        if read.is_supplementary:
            if re.match('^([0-9]+H)', read.cigarstring):
                break_count["breakcount"][read.get_blocks()[0][0]] += 1
            else:
                if len(read.get_blocks()) == 1:
                    break_count["breakcount"][read.get_blocks()[0][1] - 1] += 1
                else:
                    break_count["breakcount"][read.get_blocks()[-1][1] - 1] += 1

        if read.get_cigar_stats()[0][4] > 0:
            if re.match('^([0-9]+S)', read.cigarstring):
                break_count["breakcount"][read.get_blocks()[0][0]] += 1
            if (read.cigarstring).endswith('S'):
                if len(read.get_blocks()) == 1:
                    break_count["breakcount"][read.get_blocks()[0][1] - 1] += 1
                else:
                    break_count["breakcount"][read.get_blocks()[-1][1] - 1] += 1
    data = pd.DataFrame(break_count)
    data['position'] = data.index + 1
    data['contig'] = ref
    data = data.loc[data['breakcount'] > 0, ]
    return data


def window_break_cal(data):
    data['start_pos'] = [math.floor(x) * 100 + 300 for x in (data['position'] - 300) / 100]
    data = data.loc[data['start_pos'] >= 300, ]
    data['read_breakpoint_ratio'] = data['read_breakpoint_count'] / \
        data['read_count']
    data['index'] = data['contig'] + '_' + \
        [str(int(x)) for x in data['start_pos']]
    grouped = data.groupby(['index'])
    read_break_ratio = pd.DataFrame(grouped['read_breakpoint_ratio'].max())
    read_break_ratio['contig'] = ['_'.join(x.split("_")[:-1]) for x in read_break_ratio.index]
    read_break_ratio['start_pos'] = [int(x.split("_")[-1]) for x in read_break_ratio.index]
    read_break_ratio.index = range(read_break_ratio.shape[0])
    return read_break_ratio


def read_breakpoint_cal(args):
    if os.path.exists(os.path.join(args.output,
                                   "temp/read_breakpoint/read_breakpoint_per_window.txt")):
        return 0

    if os.path.exists(os.path.join(args.output,
                                   "temp/read_breakpoint/read_breakpoint_per_base.txt")):
        read_breakpoint_data = pd.read_csv(os.path.join(args.output,
                                                        "temp/read_breakpoint/read_breakpoint_per_base.txt"), sep="\t", index_col=0)
        window_read_breakpoint_data = window_break_cal(read_breakpoint_data)
        window_read_breakpoint_data.to_csv(os.path.join(args.output,
                                                        "temp/read_breakpoint/read_breakpoint_per_window.txt"),sep="\t")
        return 0

    samfile = pysam.AlignmentFile(args.bamfile, "rb")
    references = samfile.references
    lengths = samfile.lengths
    read_breakpoint_pool = {"contig": [],
                            "position": [],
                            "read_breakpoint_count": [],
                            "read_count": []}

    for ref, lens in zip(references, lengths):
        if lens < args.min_length:
            continue
        contig_break_data = read_breakpoint_per_contig(samfile, ref, lens)
        if contig_break_data.shape[0] > 0:
            read_breakpoint_pool["read_breakpoint_count"].extend(
                list(contig_break_data['breakcount']))
            read_breakpoint_pool["read_count"].extend(
                list(contig_break_data['readcount']))
            read_breakpoint_pool["contig"].extend(
                [ref] * contig_break_data.shape[0])
            read_breakpoint_pool["position"].extend(
                list(contig_break_data['position']))
    read_breakpoint_data = pd.DataFrame(read_breakpoint_pool)
    read_breakpoint_data.to_csv(os.path.join(args.output,
                                             "temp/read_breakpoint/read_breakpoint_per_base.txt"), sep="\t")
    window_read_breakpoint_data = window_break_cal(read_breakpoint_data)
    window_read_breakpoint_data.to_csv(os.path.join(args.output,
                                                    "temp/read_breakpoint/read_breakpoint_per_window.txt"), sep="\t")


def pileupfile_parse(args):
    """
    process pileup file
    """
    if os.path.exists(os.path.join(args.output,
                                   "temp/pileup/pileup_feature.txt")):
        return 0
    if not os.path.exists(args.pileup):
        if os.path.exists(os.path.join(args.output,
                                       "temp/pileup/contigs_pipelup.out")):
            args.pileup = os.path.join(args.output, "temp/pileup/contigs_pipelup.out")
        else:
            if not os.path.exists(args.assemblies):
                if os.path.exists(os.path.join(args.output,
                                               "temp/contig/filtered_contigs.fa")):
                    args.assemblies = os.path.join(args.output, "temp/contig/filtered_contigs.fa")
                else:
                    sys.stderr.write(f"Error: Can not find assemblies:{args.assemblies}!\n")
                    sys.exit(1)

            os.makedirs(os.path.join(args.output, "temp/pileup"), exist_ok=True)
            pileup_command = ' '.join([args.samtools,
                                       'mpileup -C 50 -A -f',
                                       args.assemblies,
                                       args.bamfile,
                                       " | awk", "'", "$3 !=", "\"N\"", "'", ">",
                                       os.path.join(args.output, "temp/pileup/contigs_pipelup.out")])
            args.pileup = os.path.join(args.output, "temp/pileup/contigs_pipelup.out")
            os.system(pileup_command)
    samfile = pysam.AlignmentFile(args.bamfile, "rb")
    contig_len = contig_pool(samfile)

    prev_contig = None
    pileup_dict = {"contig": [], "correct": [], "ambiguous": [], "insert": [],
                   "deletion": [], "disagree": [], "depth": []}
    window_pileup_dict = {"contig": [], "start_pos": [], "correct_portion": [], "ambiguous_portion": [], "disagree_portion": [],
                          "deletion_portion": [], "insert_portion": [], "normalized_coverage": [], "normalized_deviation": [], "mean_coverage": []}

    for line in open(args.pileup, "r"):
        record = line.strip().split('\t')
        if contig_len[record[0]] < args.min_length:
            continue
        if prev_contig is None:
            prev_contig = record[0]
        if record[0] != prev_contig:
            window_data = pileup_window_cal(pileup_dict)
            mean_cov = np.mean(window_data["coverage"])
            window_pileup_dict["contig"].extend(window_data["contig"])
            window_pileup_dict["start_pos"].extend(window_data["start_pos"])
            window_pileup_dict["correct_portion"].extend(
                window_data["correct_portion"])
            window_pileup_dict["ambiguous_portion"].extend(
                window_data["ambiguous_portion"])
            window_pileup_dict["disagree_portion"].extend(
                window_data["disagree_portion"])
            window_pileup_dict["deletion_portion"].extend(
                window_data["deletion_portion"])
            window_pileup_dict["insert_portion"].extend(
                window_data["insert_portion"])
            window_pileup_dict["normalized_coverage"].extend(
                window_data["coverage"] / mean_cov)
            window_pileup_dict["normalized_deviation"].extend(
                window_data["deviation"])
            window_pileup_dict["mean_coverage"].extend(
                [mean_cov] * len(window_data["start_pos"]))
            pileup_dict = {"contig": [],
                           "correct": [],
                           "ambiguous": [],
                           "insert": [],
                           "deletion": [],
                           "disagree": [],
                           "depth": []}
            prev_contig = record[0]
        pileup_dict['contig'].append(record[0])
        match_detail = record[4]
        pileup_dict['correct'].append(match_detail.count('.') + match_detail.count(','))
        pileup_dict['ambiguous'].append(match_detail.count('*'))
        pileup_dict['insert'].append(match_detail.count("+"))
        pileup_dict['deletion'].append(match_detail.count("-"))
        pileup_dict['depth'].append(int(record[3]))
        st = ''.join(re.split(r'[\+|\-][0-9]+[ATCGatcg]+', match_detail))
        numd = st.count('a') + st.count('A') + st.count('t') + st.count('T') + \
            st.count('c') + st.count('C') + st.count('g') + st.count('G')
        pileup_dict['disagree'].append(numd)

    if not os.path.exists(os.path.join(args.output, "temp/pileup")):
        os.makedirs(os.path.join(args.output, "temp/pileup"), exist_ok=True)

    data = pd.DataFrame(window_pileup_dict)
    data.to_csv(os.path.join(args.output, "temp/pileup/pileup_feature.txt"), sep="\t")
    return data


def read_cal(args, mu, dev):
    if os.path.exists(os.path.join(args.output,
                                   "temp/read_feature/read_feature.txt")):
        return 0
    samfile = pysam.AlignmentFile(args.bamfile, "rb")
    references = samfile.references
    lengths = samfile.lengths
    read_dicts = {"contig": [], "start_pos": [], "read_count": [], "proper_read_count": [], "inversion_read_count": [],
                  "clipped_read_count": [], "supplementary_read_count": [], "discordant_size_count": [], "discordant_loc_count": [], "length": []}
    for ref, lens in zip(references, lengths):
        if lens < args.min_length:
            continue
        contig_reads = samfile.fetch(ref)
        read_dict = window_read_cal(contig_reads, mu, dev)
        read_dicts["start_pos"].extend(read_dict["start_pos"])
        read_dicts["contig"].extend([ref] * len(read_dict["start_pos"]))
        read_dicts["read_count"].extend(read_dict["read_count"])
        read_dicts["proper_read_count"].extend(read_dict["proper_read_count"])
        read_dicts["inversion_read_count"].extend(
            read_dict["inversion_read_count"])
        read_dicts["clipped_read_count"].extend(
            read_dict["clipped_read_count"])
        read_dicts["supplementary_read_count"].extend(
            read_dict["supplementary_read_count"])
        read_dicts["discordant_size_count"].extend(
            read_dict["discordant_size_count"])
        read_dicts["discordant_loc_count"].extend(
            read_dict["discordant_loc_count"])
        read_dicts["length"].extend([lens] * len(read_dict["start_pos"]))
    data = pd.DataFrame(read_dicts)
    data.to_csv(os.path.join(args.output,
                             "temp/read_feature/read_feature.txt"), sep="\t")


def fragment_cal(args, mu, dev):
    if os.path.exists(os.path.join(args.output,
                                   "temp/coverage/fragment_coverage.txt")):
        return 0
    samfile = pysam.AlignmentFile(args.bamfile, "rb")
    references = samfile.references
    lengths = samfile.lengths
    frag_dict = {
        "contig": [],
        "start_pos": [],
        "normalized_fragment_coverage": [],
        "normalized_fragment_deviation": []}
    for ref, lens in zip(references, lengths):
        if lens < args.min_length:
            continue
        reads = samfile.fetch(ref)
        frag_coverage = fragment_coverage_cal(reads, mu, dev, lens)
        fragcov = window_frag_cal(frag_coverage)
        frag_dict["contig"].extend([ref] * len(fragcov['pos']))
        frag_dict["start_pos"].extend(fragcov["pos"])
        frag_dict["normalized_fragment_coverage"].extend(
            fragcov["coverage"] / np.mean(fragcov["coverage"]))
        frag_dict["normalized_fragment_deviation"].extend(fragcov["deviation"])
    data = pd.DataFrame(frag_dict)
    data.to_csv(os.path.join(args.output,
                             "temp/coverage/fragment_coverage.txt"), sep="\t")


def KAD_cal(args):
    if os.path.exists(os.path.join(args.output,
                                   "temp/KAD/KAD_window_data.txt")):
        return 0

    contig_data = pd.read_csv(os.path.join(args.output,
                                           "temp/split/contig_name.txt"), header=None)
    split_data = pd.read_csv(os.path.join(args.output,
                                          "temp/split/split_file_name.txt"), header=None)
    data = pd.concat([contig_data, split_data], axis=1)
    data.columns = ['contig', 'file']
    data.index = data['contig']
    contig_file = data.loc[:, 'file'].to_dict()

    os.makedirs(os.path.join(args.output, 'temp/KAD/temp'), exist_ok=True)
    os.makedirs(os.path.join(args.output, 'temp/KAD/KAD_data'), exist_ok=True)

    pool = multiprocessing.Pool(processes=args.threads)
    samfile = pysam.AlignmentFile(args.bamfile, "rb")
    contig_len = contig_pool(samfile)
    for contig, file in contig_file.items():
        if contig_len[contig] < args.min_length:
            continue
        try:
            t = pool.apply_async(func=KAD, args=(args, contig, file,))
        except BaseException:
            continue
    pool.close()
    pool.join()
    KAD_dict = KAD_feature(args)
    KAD_window_data = pd.DataFrame(KAD_dict)
    KAD_window_data.to_csv(os.path.join(args.output,
                                        "temp/KAD/KAD_window_data.txt"), sep="\t")

def extract_features(args):
    os.makedirs(os.path.join(args.output, 'temp', 'read_feature'), exist_ok=True)
    os.makedirs(os.path.join(args.output, 'temp', 'coverage'), exist_ok=True)
    os.makedirs(os.path.join(args.output, 'temp', 'pileup'), exist_ok=True)
    os.makedirs(os.path.join(args.output, 'temp', 'read_breakpoint'), exist_ok=True)

    samfile = pysam.AlignmentFile(args.bamfile, "rb")
    size_freq = fragment_distribution(samfile)
    mu, dev = FragMAD(size_freq)
    pool = [multiprocessing.Process(target=read_cal, args=(args, mu, dev,)),
            multiprocessing.Process(
        target=fragment_cal, args=(
            args, mu, dev,)),
            multiprocessing.Process(target=pileupfile_parse, args=(args,)),
            multiprocessing.Process(target=read_breakpoint_cal, args=(args,)),
            multiprocessing.Process(target=split_sam, args=(args,))]
    for t in pool:
        t.start()
    for t in pool:
        t.join()
    KAD_cal(args)


