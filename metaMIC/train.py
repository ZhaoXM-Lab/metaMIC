#!/usr/bin/env python

import pandas as pd
import numpy as np
import multiprocessing
import argparse,operator,os,random,sys,time
import random,subprocess
import warnings
import joblib
from sklearn.ensemble import RandomForestClassifier

base_path = os.path.split(__file__)[0]

def parseargs():
    parser=argparse.ArgumentParser(description="Calculate read features")
    parser.add_argument('--data',help='trainning data file')
    parser.add_argument('--label',help='label for train data')
    parser.add_argument('--train',help='name of training models')
    parser.add_argument("--thread",type=int, default=8,help='length of k-mer')
    args=parser.parse_args()
    return args
    
def trainforest(data,path,n):
    positive = data.loc[data['label']==1,]
    negative = data.loc[data['label']==0,]
    sub_sample = negative.iloc[random.sample(range(negative.shape[0]),positive.shape[0]),]
    balanced = pd.concat([positive,sub_sample])
    train_data = balanced.loc[:,balanced.columns != 'label']
    train_label = np.array([1 if i == 1 else -1 for i in balanced['label']]).reshape(balanced.shape[0], 1)
    rf = RandomForestClassifier(n_estimators=2000,class_weight='balanced')
    rf.fit(train_data,train_label)
    joblib.dump(rf, os.path.join(path, 'RF{}.pkl'.format(str(n))))

def training(args,data):
    model_path = os.path.join(base_path,"model",args.assembler)
    os.makedirs(model_path, exist_ok=True)
    pool=multiprocessing.Pool(processes=args.threads)
    for i in range(10):
        t = pool.apply_async(func=trainforest,args=(data,model_path,i,))
    pool.close()
    pool.join()

def train(data, args):
    train_data = pd.read_csv(data,sep="\t",index_col=0)
    train_label = pd.read_csv(args.label,sep="\t",header=None,index_col=0)
    train_data['label'] = list(train_label.loc[train_data.index,1])
    features = ['coverage_width', 'deviation_width', 'normalized_deviation','window_cov_dev', 'fragment_width', 'fragment_deviation_width','normalized_fragment_deviation', 'window_frag_cov_dev','proper_read_ratio', 'clipped_read_ratio', 'supplementary_read_ratio','inversion_read_ratio', 'discordant_loc_ratio', 'discordant_size_ratio','read_breakpoint_ratio', 'proper_read_width', 'clipped_read_width','supplementary_read_width', 'inversion_read_width','discordant_loc_width', 'discordant_size_width', 'read_breakpoint_max','disagree_width', 'correct_portion', 'ambiguous_portion','insert_portion', 'deletion_portion', 'disagree_portion', 'mean_KAD','abnormal_KAD_ratio', 'dev_KAD', 'KAD_width', 'coverage_diff','label']
    train_data = train_data.loc[:,features]
    training(args,train_data)



