# coding=utf-8
import numpy as np
import pandas as pd

seq2num_dict = dict(zip(list('ATCG'), ['01','10', '00', '11']))
seq2num_dict.update(dict(zip(list('atcg'), ['01','10', '00', '11'])))
"假设序列中不存在N/n"
"""
1. 序列清洗，去除两端包含N的序列
2. 序列打断：包含N的地方，序列需打段，去掉打断后序列低于30的read
4. 获得所有输入序列
3. cd-hit对序列进行聚类
"""


def seq2num(seq):
    new_seq = ''.join(seq2num_dict[x] for x in seq)
    return np.array([x for x in new_seq], dtype=int)

def join_two_seqs(seq1, seq2, overlap_percent=0.8, mismatch_num=2):
    """
    :param seq1: 碱基序列，必须不长于seq2
    :param seq2: 碱基序列，长度必须大于等于seq1
    :param overlap_percent: 以seq1为准
    :return:
    """
    seq1 = seq2num(seq1)
    seq2 = seq2num(seq2)
    mismatch_num = mismatch_num
    if (seq1 - seq2).__abs__().sum() <= 1:

        return seq1
    link_site = int(len(seq1)/2*overlap_percent)
    if (seq1[link_site:] - seq2[:-link_site]).sum() == 0:
        return np.concatenate((seq1[:-link_site], seq2[-link_site:]))
    elif (seq2[link_site:] - seq1[:-link_site]).sum() == 0:
        return np.concatenate((seq2[:-link_site], seq1[-link_site:]))
    else:
        return ''

def break_read(seq):
    pass
