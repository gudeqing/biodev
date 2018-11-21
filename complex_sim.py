# coding=utf-8
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from pyflow import WorkflowRunner
import glob
import time
import logging
import random
import linecache
import re
import argparse
from simulate import batch_simulate, run_pipeline, wait_bam_born, filter_vcf, sampling
__author__ = 'gdq'

"""
1. downsample RD1802566PLA 到3000x, 得到sample， S1
2. 对S1进行pipeline分析，顺便得到S1-bam
3. 基于S1-bam进行突变模拟, 得到 sample， SM
4. 重复2-3，一次模拟不同的突变类型
"""


def run_simulate(sampled_dir, batch_id, panel, sample_type, pipeline_out, vcf_file, mu_type, af_value, pyflow=True, skip_pipeline=False):
    """run pipeline to get bam, and simulate based on resulted bam 仅支持模拟一个样本"""
    if not skip_pipeline:
        run_pipeline(sampled_dir, batch_id=batch_id, panel=panel, sample_type=sample_type, out_dir=pipeline_out)
        time.sleep(3666)  # wait for run_pipeline
        sample_id_list = [os.path.basename(x).split('_')[0] for x in glob.glob(sampled_dir + '/*_R1_*.fastq*')]
        bam_dir_list = (pipeline_out + '/{}/bam/'.format(x) for x in sample_id_list )
        map(wait_bam_born, bam_dir_list)  # wait for bam.bai
    # 注： 只有当vcf仅有5列时，下面的af-value才会起到作用
    out_dir_set = batch_simulate(sampled_dir, pipeline_out=pipeline_out, vcf_list=[vcf_file],
                                 mu_type=mu_type, af_list=[af_value], pyflow=pyflow)
    return glob.glob(out_dir_set.pop() + '/*_out_reads/')[0]


def run_simulate_pipeline(sample_name, raw_dir, qc_result, depth=3000, seed=1,
                          batch_id='LK291', panel='Oncoscreen.520.v3', sample_type='FFPE',
                          snv_vcf=None, indel_vcf=None, fusion_vcf=None, cnv_vcf=None,
                          pyflow=False, run_final_pipeline=True):
    """
    累加模拟各个突变类型,仅支持模拟一个样本
    :param sample_name: 样本名
    :param raw_dir: fastq所在目录
    :param qc_result: 质控结果文件，是downsample需要的输入
    :param depth: downsample的输入，指定覆盖度
    :param seed: downsample的输入，随机取值的种子
    :param batch_id: config 信息，跑pipeline需要
    :param panel: config信息，跑pipeline需要
    :param sample_type: 样本类型，跑pipeline需要
    :param snv_vcf: 要模拟的snv信息
    :param indel_vcf: 要模拟的indel信息
    :param fusion_vcf: 要模拟的fusion信息
    :param cnv_vcf: 要模拟的cnv信息
    :param pyflow: 是否启用pyflow来进行突变模拟，若为false，则用multiprocess来管理模拟任务
    :param run_final_pipeline
    :return:
    """
    # step1: down sample
    fastq_dir = sampling(sample_name, raw_dir, qc_result, depth=depth, seeds=[seed])

    # step2: get bam and simulating based on bam
    if snv_vcf:
        fastq_dir = run_simulate(fastq_dir, batch_id, panel, sample_type, 'run_pipeline_for_down_sample',
                                 snv_vcf, "snv", 0.1, pyflow=pyflow, skip_pipeline=False)
    if indel_vcf:
        fastq_dir = run_simulate(fastq_dir, batch_id, panel, sample_type, 'run_pipeline_for_snv_simulation',
                                 indel_vcf, 'indel', 0.1, pyflow=pyflow, skip_pipeline=False)
    if fusion_vcf:
        fastq_dir = run_simulate(fastq_dir, batch_id, panel, sample_type, 'run_pipeline_for_indel_simulation',
                                 fusion_vcf, 'fusion', 0.1, pyflow=pyflow, skip_pipeline=False)
    if cnv_vcf:
        fastq_dir = run_simulate(fastq_dir, batch_id, panel, sample_type, 'run_pipeline_for_fusion_simulation',
                                 cnv_vcf, 'cnv', 0.1, pyflow=pyflow, skip_pipeline=False)
    if run_final_pipeline:
        run_pipeline(fastq_dir, batch_id=batch_id, panel=panel, sample_type=sample_type,
                     out_dir="run_pipeline_for_last_simulation")


def main():
    sample_name = 'RD1802566PLA'
    raw_dir = '/share/home/deqing.gu/simulation/tissue/'
    qc_result = '/share/home/deqing.gu/simulation/tissue/RD1802566PLA.QC.xls'
    batch_id = 'LK291'; panel = 'Oncoscreen.520.v3'; sample_type = 'FFPE'
    # exclude_vcfs = ['/share/home/deqing.gu/simulation/mutation_pos/NA12878.vcf', ]
    # target_vcfs = '/share/home/deqing.gu/simulation/mutation_pos/indel_annotation/target.indel.vcf'
    snv_vcf = '/share/home/deqing.gu/workspace/complex_simulate/simulate_info/snv.vcf'
    indel_vcf = '/share/home/deqing.gu/workspace/complex_simulate/simulate_info/indel.vcf'
    fusion_vcf = '/share/home/deqing.gu/workspace/complex_simulate/simulate_info/fusion.info'
    cnv_vcf = '/share/home/deqing.gu/workspace/complex_simulate/simulate_info/cnv.bed'
    run_simulate_pipeline(sample_name, raw_dir, qc_result, depth=3000, seed=1,
                          batch_id=batch_id, panel=panel, sample_type=sample_type,
                          snv_vcf=snv_vcf, indel_vcf=indel_vcf, fusion_vcf=fusion_vcf, cnv_vcf=cnv_vcf,
                          pyflow=True)

def main2():
    sample_name = 'RD1804206FFP'
    raw_dir = '/share/home/devdata/deqing.gu/simulate_rawdata_for_FFPE/FASTQ'
    qc_result = '/share/home/devdata/deqing.gu/simulate_rawdata_for_FFPE/QC/RD1804206FFP.QC.xls'
    batch_id = 'LK291'
    panel = 'Oncoscreen.520.v3'
    sample_type = 'FFPE'
    snv_vcf = '/share/home/devdata/deqing.gu/simulate_rawdata_for_FFPE/simulate_info/snv.vcf'
    indel_vcf = '/share/home/devdata/deqing.gu/simulate_rawdata_for_FFPE/simulate_info/indel.vcf'
    fusion_vcf = '/share/home/devdata/deqing.gu/simulate_rawdata_for_FFPE/simulate_info/fusion.info'
    cnv_vcf = '/share/home/devdata/deqing.gu/simulate_rawdata_for_FFPE/simulate_info/cnv.bed'
    run_simulate_pipeline(sample_name, raw_dir, qc_result, depth=3000, seed=1,
                          batch_id=batch_id, panel=panel, sample_type=sample_type,
                          snv_vcf=snv_vcf, indel_vcf=indel_vcf, fusion_vcf=fusion_vcf, cnv_vcf=cnv_vcf,
                          pyflow=True)

if __name__ == '__main__':
    # main()
    main2()

