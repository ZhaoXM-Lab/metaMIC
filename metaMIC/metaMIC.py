#!/usr/bin/env python

import pandas as pd
import numpy as np
import multiprocessing
from optparse import OptionParser
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
import re
import math
from Bio import SeqIO
import joblib
import warnings
from sklearn.ensemble import IsolationForest
from scipy import stats
import logging
import hashlib
import requests
import shutil
import tarfile
import gzip
from .extract import extract_features

base_path = os.path.split(__file__)[0]
contig_features = [
    'coverage_width',
    'deviation_width',
    'normalized_deviation',
    'window_cov_dev',
    'fragment_width',
    'fragment_deviation_width',
    'normalized_fragment_deviation',
    'window_frag_cov_dev',
    'proper_read_ratio',
    'clipped_read_ratio',
    'supplementary_read_ratio',
    'inversion_read_ratio',
    'discordant_loc_ratio',
    'discordant_size_ratio',
    'read_breakpoint_ratio',
    'proper_read_width',
    'clipped_read_width',
    'supplementary_read_width',
    'inversion_read_width',
    'discordant_loc_width',
    'discordant_size_width',
    'read_breakpoint_max',
    'disagree_width',
    'correct_portion',
    'ambiguous_portion',
    'insert_portion',
    'deletion_portion',
    'disagree_portion',
    'mean_KAD',
    'abnormal_KAD_ratio',
    'dev_KAD',
    'KAD_width',
    'coverage_diff',
    'length']
window_features = [
    'correct_portion',
    'ambiguous_portion',
    'disagree_portion',
    'deletion_portion',
    'insert_portion',
    'mean_KAD',
    'abnormal_KAD_ratio',
    'dev_KAD',
    'normalized_fragment_coverage',
    'normalized_fragment_deviation',
    'normalized_coverage',
    'normalized_deviation',
    'mean_coverage',
    'coverage_diff',
    'read_breakpoint_ratio',
    'proper_read_ratio',
    'inversion_read_ratio',
    'clipped_read_ratio',
    'supplementary_read_ratio',
    'discordant_loc_ratio',
    'discordant_size_ratio']


def get_opts(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=' Reference-free Misassembly Identification and Correction of metagenomic assemblies')

    subparsers = parser.add_subparsers(title='metaMIC subcommands',
                                       dest='cmd',
                                       metavar='')

    download_model = subparsers.add_parser('download_model',
                                           help='download model')

    extract_feature = subparsers.add_parser('extract_feature',
                                            help='Extract features from inputs.')

    extract_feature.add_argument(
        "-t",
        "--threads",
        dest="threads",
        required=False,
        type=int,
        default=8,
        help='Maximum number of threads [default: 8]')

    extract_feature.add_argument(
        "--bam",
        dest="bamfile",
        required=False,
        help="index bam file for alignment")

    extract_feature.add_argument(
        "--r1",
        dest="read1",
        required=False,
        help="read1")

    extract_feature.add_argument(
        "--r2",
        dest="read2",
        required=False,
        help="read2")

    extract_feature.add_argument(
        "-p",
        "--r",
        dest="read",
        help="smart pairing (ignoring #2 fasta/q)")

    extract_feature.add_argument(
        "-c",
        "--contig",
        dest="assemblies",
        required=True,
        help="fasta file of assembled contigs")

    extract_feature.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help='output directory for AIM results')

    extract_feature.add_argument(
        "--pileup",
        dest="pileup",
        default='pileup',
        required=True,
        help="path to pileup file [samtools mpileup]")

    extract_feature.add_argument(
        "-m",
        "--mode",
        dest="mode",
        required=True,
        default='meta',
        help="Applied to single genomic/metagenomic assemblies [meta/single]")

    extract_feature.add_argument(
        "-l",
        "--mlen",
        dest="min_length",
        type=int,
        default=5000,
        required=False,
        help='Minimum contig length [default: 5000bp]')

    extract_feature.add_argument(
        "--samtools",
        dest="samtools",
        default='samtools',
        required=False,
        help="path to samtools")

    extract_feature.add_argument(
        "--jellyfish",
        dest="jellyfish",
        required=False,
        default="jellyfish",
        help="path to jellyfish")



    predict = subparsers.add_parser('predict',
                                    help='Predict.')
    predict.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help='output directory for AIM results')

    predict.add_argument(
        "-m",
        "--mode",
        dest="mode",
        required=True,
        default='meta',
        help="Applied to single genomic/metagenomic assemblies [meta/single]")

    predict.add_argument(
        "-c",
        "--contig",
        dest="assemblies",
        required=True,
        help="fasta file of assembled contigs")

    predict.add_argument(
        "-a",
        "--assembler",
        dest="assembler",
        required=False,
        default="MEGAHIT",
        help="The assembler-specific model or user-trained model used for assembled fasta file [MEGAHIT/IDBA_UD/[new training model specified by users]]")

    predict.add_argument(
        "--st",
        dest="score_thred",
        required=False,
        type=float,
        default=None,
        help="Threshold of metaMIC contig score for correcting misassemblies in metagenomics[default:0.8]"
    )

    predict.add_argument(
        "-l",
        "--mlen",
        dest="min_length",
        type=int,
        required=False,
        default=5000,
        help='Minimum contig length [default: 5000bp]')

    predict.add_argument(
        "-s",
        "--slen",
        dest="split_length",
        type=int,
        required=False,
        default=1000,
        help='Minimum length of splitted fragments [default: 1000bp]')

    predict.add_argument(
        "--nb",
        dest="break_count",
        type=int,
        required=False,
        default=5,
        help='Threshold of read breakpoint counts for correcting misassemblies in metagenomics')

    predict.add_argument(
        "--rb",
        dest="break_ratio",
        type=float,
        required=False,
        default=0.1,
        help='Threshold of read breakpoint ratio for correcting misassemblies in metagenomics')

    predict.add_argument(
        "--at",
        dest="anomaly_thred",
        required=False,
        type=float,
        default=0.9,
        help='Threshold of anomaly score for correcting misassemblies in metagenomics')


    train = subparsers.add_parser('train',
                                  help='Train model.')

    train.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help='output directory for AIM results')

    train.add_argument(
        "--label",
        dest="label",
        help='Misassembly label of contigs for training assemblies')

    train.add_argument(
        "-a",
        "--assembler",
        dest="assembler",
        default="MEGAHIT",
        help="The name of the directory of the trained model.")

    train.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=8,
        help='Maximum number of CPUs [default: 8]')



    if not args:
        parser.print_help(sys.stderr)
        sys.exit()

    return parser.parse_args(args)

def get_file_md5(fname):
    """
    Calculate Md5 for downloaded file
    """
    m = hashlib.md5()
    with open(fname,'rb') as fobj:
        while True:
            data = fobj.read(4096)
            if not data:
                break
            m.update(data)

    return m.hexdigest()

def download_model(path, MD5):
    os.makedirs(os.path.join(base_path, 'model'),exist_ok=True)
    download_path = os.path.join(base_path, 'model',os.path.split(path)[1])
    print(download_path)
    with requests.get(path, stream=True) as r:
        with open(download_path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)
    print('Download finished. Checking MD5...')
    if os.path.exists(os.path.join(base_path, 'model', 'MEGAHIT')):
        shutil.rmtree(os.path.join(base_path, 'model', 'MEGAHIT'))

    if os.path.exists(os.path.join(base_path, 'model', 'IDBA_UD')):
        shutil.rmtree(os.path.join(base_path, 'model', 'IDBA_UD'))

    if os.path.exists(os.path.join(base_path, 'model', 'metaSPAdes')):
        shutil.rmtree(os.path.join(base_path, 'model', 'metaSPAdes'))


    if get_file_md5(download_path) == MD5:
        try:
            tar = tarfile.open(download_path, "r:gz")
            file_names = tar.getnames()
            for file_name in file_names:
                tar.extract(file_name, os.path.join(base_path, 'model'))
            tar.close()
        except Exception:
            sys.stderr.write(
                f"Error: cannot unzip the file.")
            sys.exit(1)
        os.remove(download_path)

        model_file = os.listdir(os.path.join(base_path, 'model', os.path.split(path)[1].split('.')[0]))
        for temp_gz in model_file:
            try:
                g_file = gzip.GzipFile(os.path.join(os.path.join(base_path, 'model', os.path.split(path)[1].split('.')[0]),temp_gz))
                open(os.path.join(os.path.join(base_path, 'model', os.path.split(path)[1].split('.')[0]),temp_gz[:-3]), "wb+").write(g_file.read())
                g_file.close()
            except Exception:
                sys.stderr.write(
                    f"Error: cannot unzip the file.")
                sys.exit(1)
            os.remove(os.path.join(os.path.join(base_path, 'model', os.path.split(path)[1].split('.')[0]),temp_gz))

    else:
        os.remove(download_path)
        sys.stderr.write(
            f"Error: MD5 check failed, removing '{download_path}'.\n")
        sys.exit(1)

def download():
    print('Downloading model for MEGAHIT')
    download_model('https://zenodo.org/record/4781819/files/MEGAHIT.tar.gz', 'da9038af3582eea04288775a72003e6b')
    print('Downloading model for IDBA_UD')
    download_model('https://zenodo.org/record/4781819/files/IDBA_UD.tar.gz', '9f6835d3033a177055343ecaa78889bc')
    print('Downloading model for metaSPAdes.')
    download_model('https://zenodo.org/record/5768805/files/metaSPAdes.tar.gz',
                   'md5:0c759a28d0ba3490eeba098395c88586 ')

def filter_contig(options):
    """
    remove assemblies with length smaller than required minimal length
    """
    input = SeqIO.parse(options.assemblies, "fasta")
    input = (record for record in input if len(record.seq) > 1000)
    os.makedirs(os.path.join(options.output, "temp/contig"), exist_ok=True)
    SeqIO.write(input, os.path.join(options.output, "temp/contig/filtered_contigs.fa"), "fasta")
    options.assemblies = os.path.join(options.output, "temp/contig/filtered_contigs.fa")
    return options


def mapping(options):
    """
    mapping reads to assemblies
    """
    command_line = options.bwa + ' index ' + options.assemblies
    os.system(command_line)
    os.makedirs(os.path.join(options.output, "temp/sam/"), exist_ok=True)
    output_bam = os.path.join(options.output, "temp/sam/contigs.filter.sort.bam")
    options.bamfile = output_bam
    if options.read is None:
        command_bwa = options.bwa + ' mem -a -t ' + str(options.threads) + ' ' + options.assemblies + ' ' + options.read1 + ' ' + \
            options.read2 + ' | ' + options.samtools + ' view -h -q 10 -m 50 -F 4 -b | ' + \
            options.samtools + ' sort > ' + output_bam
    else:
        command_bwa = options.bwa + ' mem -a -p -t ' + str(options.threads) + ' ' + options.assemblies + ' ' + options.read + \
            ' | ' + options.samtools + ' view -h -q 10 -m 50 -F 4 -b | ' + \
            options.samtools + ' sort > ' + output_bam
    os.system(command_bwa)
    command_index = options.samtools + ' index ' + output_bam
    os.system(command_index)
    return options


def extract_feature(options):
    """
    Extract four types of features from bamfiles and generate contig-based/window-based feature matrix
    """
    # check if feature file exists
    status = feature_exist(options)
    if status:
        print("feature files exist and will not re-extract features")
    else:
        extract_features(options)
        # check output features
        check_feature(options)
    # Generate contig-based/window-based matrix
    if options.mode == 'single':
        window_matrix = cal_feature(options)
    else:
        window_matrix, contig_matrix = cal_feature(options)


def cov_thread_cal(data):
    up = np.percentile(np.array(data), 95)
    low = np.percentile(np.array(data), 5)
    return up, low


def contig_fea_generate(data):
    # Generate contig-based feature matrix
    grouped = data.groupby(['contig'])
    mean_data = grouped.mean()
    mean_data = mean_data.loc[:, mean_data.columns != 'start_pos']
    data['count'] = 1
    data['proper_status'] = data['proper_read_ratio'] <= 0.9
    data['clipped_status'] = data['clipped_read_ratio'] >= 0.1
    data['supplementary_status'] = data['supplementary_read_ratio'] >= 0.1
    data['inversion_status'] = data['inversion_read_ratio'] >= 0.1
    data['discordant_loc_status'] = data['discordant_loc_ratio'] >= 0.1
    data['discordant_size_status'] = data['discordant_size_ratio'] >= 0.1
    data['KAD_status'] = data['mean_KAD'] > 0.75
    data['disagree_status'] = data['disagree_portion'] > 0.1
    data['coverage_status'] = (data['normalized_coverage'] > cov_thread_cal(list(data['normalized_coverage']))[0]) + \
                              (data['normalized_coverage'] < cov_thread_cal(list(data['normalized_coverage']))[1])
    data['fragment_status'] = (data['normalized_fragment_coverage'] > cov_thread_cal(list(data['normalized_fragment_coverage']))[0]) + \
                              (data['normalized_fragment_coverage'] < cov_thread_cal(list(data['normalized_fragment_coverage']))[1])
    data['deviation_status'] = data['normalized_deviation'] > cov_thread_cal(list(data['normalized_deviation']))[0]
    data['fragment_deviation_status'] = data['normalized_fragment_deviation'] > cov_thread_cal(
        list(data['normalized_fragment_deviation']))[0]
    grouped = data.groupby(['contig'])
    width_data = pd.concat([grouped['proper_status'].sum() /
                            grouped['count'].sum(), grouped['clipped_status'].sum() /
                            grouped['count'].sum(), grouped['supplementary_status'].sum() /
                            grouped['count'].sum(), grouped['inversion_status'].sum() /
                            grouped['count'].sum(), grouped['discordant_loc_status'].sum() /
                            grouped['count'].sum(), grouped['discordant_size_status'].sum() /
                            grouped['count'].sum(), grouped['KAD_status'].sum() /
                            grouped['count'].sum(), grouped['disagree_status'].sum() /
                            grouped['count'].sum(), grouped['coverage_status'].sum() /
                            grouped['count'].sum(), grouped['fragment_status'].sum() /
                            grouped['count'].sum(), grouped['deviation_status'].sum() /
                            grouped['count'].sum(), grouped['fragment_deviation_status'].sum() /
                            grouped['count'].sum()], axis=1)
    width_data.columns = [
        'proper_read_width',
        'clipped_read_width',
        'supplementary_read_width',
        'inversion_read_width',
        'discordant_loc_width',
        'discordant_size_width',
        'KAD_width',
        'disagree_width',
        'coverage_width',
        'fragment_width',
        'deviation_width',
        'fragment_deviation_width']
    dev_data = pd.concat([pd.DataFrame(np.sqrt(grouped['normalized_coverage'].var())),
                          pd.DataFrame(np.sqrt(grouped['normalized_fragment_coverage'].var()))], axis=1)
    dev_data.columns = ['window_cov_dev', 'window_frag_cov_dev']
    break_data = pd.DataFrame(grouped['read_breakpoint_ratio'].max())
    break_data.columns = ['read_breakpoint_max']
    contig_data = pd.concat([mean_data, width_data, dev_data, break_data], axis=1)
    contig_data = contig_data.loc[:, contig_features]
    contig_data = contig_data.fillna(0)
    contig_data = contig_data.replace(np.inf, 0)
    return contig_data


def cal_feature(options):
    read_feature = pd.read_csv(os.path.join(options.output,
                                            "temp/read_feature/read_feature.txt"), sep="\t", index_col=0)
    read_feature['proper_read_ratio'] = read_feature['proper_read_count'] / read_feature['read_count']
    read_feature['inversion_read_ratio'] = read_feature['inversion_read_count'] / read_feature['proper_read_count']
    read_feature['clipped_read_ratio'] = read_feature['clipped_read_count'] / read_feature['proper_read_count']
    read_feature['supplementary_read_ratio'] = read_feature['supplementary_read_count'] / read_feature['proper_read_count']
    read_feature['discordant_loc_ratio'] = read_feature['discordant_loc_count'] / read_feature['proper_read_count']
    read_feature['discordant_size_ratio'] = read_feature['discordant_size_count'] / read_feature['proper_read_count']
    window_read_data = read_feature.loc[:,['contig',
                                         'start_pos',
                                         'proper_read_ratio',
                                         'inversion_read_ratio',
                                         'clipped_read_ratio',
                                         'supplementary_read_ratio',
                                         'discordant_loc_ratio',
                                         'discordant_size_ratio',
                                         'length']]

    frag_coverage = pd.read_csv(os.path.join(options.output,
                                             "temp/coverage/fragment_coverage.txt"), sep="\t", index_col=0)
    pileup_feature = pd.read_csv(os.path.join(options.output,
                                              "temp/pileup/pileup_feature.txt"), sep="\t", index_col=0)
    KAD_feature = pd.read_csv(os.path.join(options.output,
                                           "temp/KAD/KAD_window_data.txt"), sep="\t", index_col=0)
    breakpoint_data = pd.read_csv(os.path.join(options.output,
                                               "temp/read_breakpoint/read_breakpoint_per_window.txt"), sep="\t", index_col=0)
    window_data = pd.merge(window_read_data, frag_coverage, on=['contig', 'start_pos'])
    window_data = pd.merge(window_data, pileup_feature, on=['contig', 'start_pos'])
    window_data = pd.merge(window_data, KAD_feature, on=['contig', 'start_pos'])
    window_data = pd.merge(window_data, breakpoint_data, on=['contig', 'start_pos'], how='left')
    window_data = window_data.fillna(0)
    window_data['coverage_diff'] = window_data['normalized_coverage'] - \
        window_data['normalized_fragment_coverage']
    os.makedirs(os.path.join(options.output, "feature_matrix"), exist_ok=True)
    window_data = window_data.loc[window_data['mean_coverage'] > 5, ]
    window_data.to_csv(os.path.join(options.output, "feature_matrix/window_fea_matrix.txt"), sep="\t")
    if options.mode == 'single':
        return window_data
    else:
        contig_data = contig_fea_generate(window_data)
        contig_data.to_csv(os.path.join(options.output,
                                        "feature_matrix/contig_fea_matrix.txt"), sep="\t")
        return window_data, contig_data


def feature_exist(options):
    status = 1

    def check_status(f):
        if not os.path.exists(f):
            nonlocal status
            status = 0

    check_status(os.path.join(options.output,
                                       "temp/read_feature/read_feature.txt"))
    check_status(os.path.join(options.output,
                                       "temp/coverage/fragment_coverage.txt"))
    check_status(os.path.join(options.output,
                                       "temp/read_breakpoint/read_breakpoint_per_window.txt"))
    check_status(os.path.join(options.output,
                                       "temp/pileup/pileup_feature.txt"))
    check_status(os.path.join(options.output,
                                       "temp/KAD/KAD_window_data.txt"))

    return status


def check_feature(options):
    def check_path(f):
        if not os.path.exists(f):
            sys.stderr.write(f"Error: Expected file '{f}' does not exist\n")
            sys.exit(1)
    check_path(os.path.join(options.output,
                            "temp/read_feature/read_feature.txt"))
    check_path(os.path.join(options.output,
                            "temp/coverage/fragment_coverage.txt"))
    check_path(os.path.join(options.output,
                            "temp/read_breakpoint/read_breakpoint_per_window.txt"))
    check_path(os.path.join(options.output,
                            "temp/pileup/pileup_feature.txt"))
    check_path(os.path.join(options.output,
                            "temp/KAD/KAD_window_data.txt"))


def predict(options, data):
    # identify misassembled metagenomic contigs
    min_length = options.min_length
    test_data = data.loc[data['length'] > min_length, data.columns != 'length']
    modelfilename = '/'.join(["model", options.assembler])
    model_path = os.path.join(base_path, modelfilename)
    if not os.path.exists(model_path):
        f = model_path
        sys.stderr.write(
            f"Error: Expected training model '{f}' does not exist\n")
        sys.exit(1)
    score = pd.DataFrame(np.zeros([test_data.shape[0], 10]))
    score.index = test_data.index
    for i in range(10):
        rf = joblib.load(model_path + "/RF" + str(i) + '.pkl')
        pro = pd.DataFrame(rf.predict_proba(test_data))
        pro.index = test_data.index
        score.loc[pro.index, i] = pro[1]
    score['metaMIC_contig_score'] = score.mean(axis=1)
    score = score.loc[:, ['metaMIC_contig_score']]
    score['length'] = data.loc[score.index, 'length']
    os.makedirs(options.output, exist_ok=True)
    score.to_csv(os.path.join(options.output, "metaMIC_contig_score.txt"), sep="\t")
    return score


def findcut(options, score):
    # score cutoff for identify candidate misassembled contigs
    if options.assembler == 'MEGAHIT':
        score_cut = 0.8
    elif options.assembler == 'IDBA_UD':
        score_cut = 0.5
    else:
        score_cut = np.percentile(score['metaMIC_contig_score'], 95)
        score_cut = round(score_cut * 10) / 10
    return score_cut


def Isolation_forest(options, window_data):
    """
    Isolation Forest
    """
    Xdata = window_data.loc[:, window_features]
    Xdata = Xdata.fillna(0)
    Xdata = Xdata.replace(np.inf, 0)
    n_samples = Xdata.shape[0]
    if options.mode == 'meta':
        iso_parameter = {"outlier_fraction": 0.001, "outlier_threshold": 0.01}
    else:
        iso_parameter = {
            "outlier_fraction": 0.0001,
            "outlier_threshold": 0.001}
    # fit the model
    score_pred = np.zeros([n_samples])
    for i in range(5):
        clf = IsolationForest(contamination=iso_parameter["outlier_fraction"])
        clf.fit(Xdata)
        score_pred += clf.decision_function(Xdata)
    score_pred = score_pred / 5
    Xdata['anomaly_score'] = 1 - score_pred
    threshold = stats.scoreatpercentile(
        1 - score_pred, 100 * (1 - iso_parameter["outlier_threshold"]))
    Xdata['anomaly_thred'] = threshold
    score_pred_data = Xdata.loc[:, ['anomaly_score',
                                    'anomaly_thred', 'read_breakpoint_ratio']]
    score_pred_data['start_pos'] = [
        int(x.split("_")[-1]) for x in score_pred_data.index]
    score_pred_data['contig'] = ['_'.join(x.split("_")[:-1])
                                 for x in score_pred_data.index]
    return score_pred_data


def breakpoint_detect(options, data):
    """
    Localize misassembly breakpoints in single genomic/metagenomic assemblies
    """
    if not os.path.exists(os.path.join(options.output,
                                       "temp/read_breakpoint/read_breakpoint_per_base.txt")):
        f = os.path.join(options.output,
                                       "temp/read_breakpoint/read_breakpoint_per_base.txt")
        sys.stderr.write(f"Error: Expected file '{f}' does not exist\n")
        sys.exit(1)

    read_breakpoint = pd.read_csv(os.path.join(options.output,
                                               "temp/read_breakpoint/read_breakpoint_per_base.txt"), sep="\t", index_col=0)
    read_breakpoint['start_pos'] = [int(x) * 100 + 300 for x in list((read_breakpoint['position'] - 300) / 100)]
    data.index = data['contig'] + "_" + [str(int(x)) for x in data['start_pos']]
    if options.mode == "single":
        score_pred_data = Isolation_forest(options, data)
        score_pred_data.to_csv(options.output + "/anomaly_score.txt", sep="\t")
        score_pred_data = score_pred_data.loc[score_pred_data['anomaly_score'] > 0.95, ]
        score_pred_data.index = range(score_pred_data.shape[0])
        result = pd.merge(
            read_breakpoint, score_pred_data, on=[
                'contig', 'start_pos'])
        result = result.iloc[np.argsort(-result['read_breakpoint_count']), ]
        result['id'] = result['contig'] + "_" + \
            [str(int(x)) for x in result['start_pos']]
        result = result.drop_duplicates(['id'])
        result = result.loc[result['read_breakpoint_count'] /
                            result['read_count'] > 0.2, ]
        result = result.loc[result['read_breakpoint_count'] > 5, ]
        contig_len = data.loc[:, ['contig', 'length']].drop_duplicates()
        contig_len.index = contig_len['contig']
        result['contig_length'] = list(
            contig_len.loc[result['contig'], 'length'])
        breakpoint_result = result.loc[:,
                                       ["contig",
                                        "position",
                                        "read_breakpoint_count",
                                        "read_breakpoint_ratio",
                                        "anomaly_score",
                                        "anomaly_thred",
                                        "contig_length"]]
        breakpoint_result.columns = [
                                    "contig",
                                    "misassembly_breakpoint",
                                    "read_breakpoint_count",
                                    "read_breakpoint_ratio",
                                    "anomaly_score",
                                    "anomaly_thred",
                                    "contig_length"]
        breakpoint_result.to_csv(os.path.join(options.output,
                                              "misassembly_breakpoint.txt"), sep="\t")
    else:
        if not os.path.exists(os.path.join(options.output, 'metaMIC_contig_score.txt')):
            f = options.contig_score
            sys.stderr.write(f"Error: Expected file '{f}' does not exist\n")
            sys.exit(1)

        contig_score = pd.read_csv(os.path.join(options.output, 'metaMIC_contig_score.txt'), sep="\t", index_col=0)
        if options.score_thred:
            score_cut = options.score_thred
        else:
            score_cut = findcut(options, contig_score)
        filtered = contig_score.loc[contig_score['metaMIC_contig_score'] > score_cut, ]
        data = data[data['contig'].isin(filtered.index)]
        score_pred_data = Isolation_forest(options, data)
        score_pred_data.index = range(score_pred_data.shape[0])
        score_pred_data.to_csv(os.path.join(options.output, "anomaly_score.txt"), sep="\t")
        score_pred_data = score_pred_data.iloc[np.argsort( -score_pred_data['anomaly_score']), ]
        score_pred_data = score_pred_data.drop_duplicates(['contig'], keep='first')
        contig_breakpoints = score_pred_data.loc[:, ['contig', 'start_pos', 'anomaly_score',
                                                     'anomaly_thred', 'read_breakpoint_ratio']]
        result = pd.merge(read_breakpoint, contig_breakpoints, on=['contig', 'start_pos'], how='right')
        result['read_breakpoint_count'] = result['read_breakpoint_count'].fillna(0)
        result['position'] = result['position'].fillna(0)
        result = result.iloc[np.argsort(-result['read_breakpoint_count']), ]
        result = result.drop_duplicates(['contig'], keep='first')
        result.loc[result['position'] == 0, 'position'] = result.loc[result['position'] == 0, 'start_pos'] + 50
        result['contig_length'] = list(filtered.loc[result['contig'], 'length'])
        result['metaMIC_contig_score'] = list(filtered.loc[result['contig'], 'metaMIC_contig_score'])
        breakpoint_result = result.loc[:,
                                       ["contig",
                                        "position",
                                        "read_breakpoint_count",
                                        "read_breakpoint_ratio",
                                        "anomaly_score",
                                        "anomaly_thred",
                                        "metaMIC_contig_score",
                                        "contig_length"]]
        breakpoint_result.columns = [
                                    'contig',
                                    'misassembly_breakpoint',
                                    'read_breakpoint_count',
                                    'read_breakpoint_ratio',
                                    'anomaly_score',
                                    "anomaly_thred",
                                    'metaMIC_contig_score',
                                    "contig_length"]
        breakpoint_result.to_csv(os.path.join(options.output,
                                              "misassembly_breakpoint.txt"), sep="\t")
    return breakpoint_result


def correct(options, breakpoint_result):
    """
    Correct misassemblies
    """
    breakpoint_result = breakpoint_result.loc[breakpoint_result['misassembly_breakpoint']
                                              > options.split_length, ]
    breakpoint_result = breakpoint_result.loc[(breakpoint_result['contig_length'] - breakpoint_result['misassembly_breakpoint'])
                                              > options.split_length, ]
    breakpoint_result = breakpoint_result.loc[breakpoint_result['read_breakpoint_count']
                                              > options.break_count, ]
    breakpoint_result = breakpoint_result.loc[breakpoint_result['read_breakpoint_ratio']
                                              > options.break_ratio, ]
    breakpoint_result = breakpoint_result.loc[breakpoint_result['anomaly_score']
                                              > options.anomaly_thred, ]
    corrected_contig_file = os.path.join(options.output, "metaMIC_corrected_contigs.fa")
    corrected_file = open(corrected_contig_file, "w")
    original_file = options.assemblies
    input = SeqIO.parse(original_file, "fasta")
    breakcontigs = list(np.unique(breakpoint_result['contig']))
    for record in input:
        if record.id in breakcontigs:
            corrected_file.write(">" + record.id + "_1\n")
            breakpoint = int(list(breakpoint_result.loc[breakpoint_result['contig'] == record.id, 'misassembly_breakpoint'])[0])
            corrected_file.write(str(record.seq[:breakpoint]) + "\n")
            corrected_file.write(">" + record.id + "_2\n")
            corrected_file.write(str(record.seq[breakpoint:]) + "\n")
        else:
            corrected_file.write(">" + record.id + "\n")
            corrected_file.write(str(record.seq) + "\n")
    corrected_file.close()
    print("A total of " + str(breakpoint_result.shape[0]) + " misassembled contigs are corrected")

from .train import train
def train_model(options, data):
    train(data, options)


def bamindex(options):
    command_index = options.samtools + ' index ' + options.bamfile
    os.system(command_index)


def validate_options(options):
    def expect_file(f):
        if f is not None:
            if not os.path.exists(f):
                sys.stderr.write(
                    f"Error: Expected file '{f}' does not exist\n")
                sys.exit(1)

    def expect_mode(f):
        if f not in ['single', 'meta']:
            options.mode = 'meta'
            logging.warning("Using default mode: [meta]")

    if options.cmd != 'download_model':
        if not os.path.exists(options.output):
            os.makedirs(options.output, exist_ok=True)

    if options.cmd == 'extract_feature':
        if os.path.exists(os.path.join(options.output, "temp/contig/filtered_contigs.fa")):
            options.assemblies = os.path.join(options.output, "temp/contig/filtered_contigs.fa")
        else:
            expect_file(options.assemblies)
            options = filter_contig(options)
        expect_file(options.pileup)
        expect_mode(options.mode)
        if options.bamfile is not None:
            expect_file(options.bamfile)
            bamindex(options)
        else:
            if options.read is None and (options.read1 is None or options.read2 is None):
                sys.stderr.write(
                    f"Error: Expected read1 and read2.\n")
                sys.exit(1)

            if options.read and (options.read1 or options.read2):
                sys.stderr.write(
                    f"Input read1/read2 or smart pairing (ignoring #2 fasta/q).\n")
                sys.exit(1)

            if os.path.exists(os.path.join(options.output, "temp/sam/contig.filter.sort.bam")):
                options.bamfile = os.path.join(options.output, "temp/sam/contigs.filter.sort.bam")
            else:
                print("Bamfile not exists, Mapping paired-end reads to assemblies")
                options = mapping(options)

        if not os.path.exists(options.bamfile):
            f = options.bamfile
            sys.stderr.write(f"Error: Expected file '{f}' does not exist\n")
            sys.exit(1)

    if options.cmd == 'train':
        if options.label is None:
            sys.stderr.write(
                "Error: misassembly labels of training contigs does not exist\n")
            sys.exit(1)
        elif not os.path.exists(options.label):
            f = options.label
            sys.stderr.write(f"Error: Expected file '{f}' does not exist\n")
            sys.exit(1)
        if options.assembler in ["MEGAHIT", "IDBA_UD"]:
            sys.stderr.write(
                "Error: Name is same as the default model, change to another name not in [MEGAHIT, IDBA_UD].\n")
            sys.exit(1)

    if options.cmd == 'predict':
        if os.path.exists(os.path.join(options.output, "temp/contig/filtered_contigs.fa")):
            options.assemblies = os.path.join(options.output, "temp/contig/filtered_contigs.fa")
        else:
            expect_file(options.assemblies)
            options = filter_contig(options)

        expect_mode(options.mode)
    return options


def main():
    args = sys.argv[1:]
    options = get_opts(args)

    warnings.filterwarnings("ignore")
    logger = logging.getLogger('metaMIC')
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
    logger.addHandler(sh)

    ########### Validate input files ###########
    logger.info('Start metaMIC')

    options = validate_options(options)
    if options.cmd == 'download_model':
        download()

    ########### Extract four types of features from bamfile ###########
    if options.cmd == 'extract_feature':
        logger.info('Step: Feature extraction')
        extract_feature(options)

    ########### Training models specified by users or not ###########
    if options.cmd == 'train':
        logger.info('Step: Training new models')

        train_model(options, os.path.join(options.output, 'feature_matrix/contig_fea_matrix.txt'))
        logger.info("Finished")
        return 0

    ########### Identify misassembled contigs [for only metagenomics] ########
    if options.cmd == 'predict':
        if options.mode == 'meta':
            logger.info('Step: Identify misassembled metagenomic contigs')
            contig_data = pd.read_csv(os.path.join(options.output, 'feature_matrix/contig_fea_matrix.txt'), sep="\t", index_col=0)
            score = predict(options, contig_data)

        logger.info('Step: Localize misassembly breakpoints')
        window_data = pd.read_csv(os.path.join(options.output, 'feature_matrix/window_fea_matrix.txt'), sep="\t", index_col=0)
        breakpoint_result = breakpoint_detect(options, window_data)
        logger.info('Step: Correct misassemblies')
        correct(options, breakpoint_result)
        logger.info("Finished")

if __name__ == '__main__':
    main()
