# coding=utf-8
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
import shlex
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns


def down_sample_bam(bam, outdir=os.getcwd(), step=5, samtools="samtools", threads=3, pool_size=4):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    percents = list(range(10, 100, step))
    outs = []
    cmds = list()
    for seed, percent in enumerate(percents):
        cmd = "{} view ".format(samtools)
        cmd += "-s {}.{} ".format(seed, percent)
        out_name = os.path.join(outdir, os.path.basename(bam)[:-3]+"{}.bam".format(percent))
        cmd += "-o {} ".format(out_name)
        cmd += "-@ {} ".format(threads)
        cmd += "{} ".format(bam)
        cmds.append(shlex.split(cmd))
        outs.append(out_name)
    target = os.path.join(outdir, os.path.basename(bam)[:-3] + "{}.bam".format(100))
    if not os.path.exists(target):
        os.symlink(bam, target)
    outs.append(target)
    with ThreadPoolExecutor(pool_size) as pool:
        pool.map(subprocess.check_call, cmds)
    return outs


def feature_count(bam_list, annotation, featureCounts="featureCounts", paired=True, threads=3, pool_size=4):
    cmds = list()
    outs = list()
    for bam in bam_list:
        cmd = "{} ".format(featureCounts)
        if paired:
            cmd += "-p "
        cmd += "-t exon -g gene_id "
        cmd += "-a {} ".format(annotation)
        cmd += "-T {} ".format(threads)
        out_name = bam + '.counts.txt'
        cmd += "-o {} ".format(out_name)
        cmd += "{} ".format(bam)
        cmds.append(shlex.split(cmd))
        outs.append(out_name)

    with ThreadPoolExecutor(pool_size) as pool:
        pool.map(subprocess.check_call, cmds)

    # for each in bam_list:
    #     os.remove(each)

    return outs


def exp_calculator_with_count(count_table_file, exp_type='both', out_prefix=None,
                              mean_read_len=149, upper_count_percent_limit=0.1):
    """
    calculate fpkm and tpm based on count table with second column containing gene length.
    :param count_table_file: example:
    -----------
    gene_id gene_length sample1 sample2
    gene1   1001    29  50
    gene2   1300    30  14
    -----------
    :param exp_type: expression type, fpkm, tpm, or 'both'. default:'both'.
    :param out_prefix: out file prefix
    :param mean_read_len: read 的平均长度
    :param upper_count_percent_limit: 基因的count占比上限，超过该上限的gene的count不计入总count
    :return: rpkm_dict, tpm_dict
    """
    if not out_prefix:
        out_prefix = count_table_file[:-4]
    if exp_type.lower() not in ['fpkm', 'tpm', 'both']:
        raise Exception('exp_type should be fpkm or tpm or both')
    count_table = pd.read_table(count_table_file, index_col=0, header=0)
    columns = count_table.columns
    gene_len = count_table[columns[0]] + 1
    if gene_len.max() > 200000:
        print("max gene length: ", gene_len.max())
        print('Warning: The maximum gene length is abnormal! But we can do nothing!')

    if gene_len.min() <= mean_read_len:
        print("min gene length: ", gene_len.min())
        sn = count_table[gene_len <= mean_read_len].shape[0]
        print("we will discard {} genes with length smaller than mean read length!".format(sn))
        count_table = count_table[gene_len > mean_read_len]
        gene_len = count_table[columns[0]] + 1

    rpkm_dict = dict()
    tpm_dict = dict()
    for sample in columns[1:]:
        # 如果某个基因的count占比超过10%，那么这个基因可能是过度扩增，计算总count时排除他
        total_counts = count_table[sample].sum()
        if total_counts <= 10:
            raise Exception("column/sample {} has nearly zero low total count".format(sample))
        count_percent = count_table[sample] / total_counts
        exclude_counts = count_table[sample][count_percent >= upper_count_percent_limit].sum()
        if exclude_counts >= 1:
            print("Exclude {} counts from total counts for normalization ".format(exclude_counts))
            total_counts = total_counts - exclude_counts

        # Divide the read counts by the length of each gene in kilobases.
        # This gives you reads per kilo-base (RPK)
        rpk = count_table[sample]/gene_len
        # get rpkm/fpkm
        if exp_type == 'fpkm' or exp_type == 'both':
            rpkm = rpk/total_counts*1000000*1000
            rpkm_dict[sample] = rpkm
        # get tpm
        if exp_type == 'tpm' or exp_type == 'both':
            norm_gene_len_total_counts = sum(rpk)
            tpm = rpk/norm_gene_len_total_counts*1000000
            tpm_dict[sample] = tpm
    # save results
    if exp_type == 'fpkm' or exp_type == 'both':
        df_rpkm = pd.DataFrame(rpkm_dict)
        df_rpkm.to_csv(out_prefix+'.fpkm.xls', sep='\t')
    if exp_type == 'tpm' or exp_type == 'both':
        df_tpm = pd.DataFrame(tpm_dict)
        df_tpm.to_csv(out_prefix+'.tpm.xls', sep='\t')
    df_count = count_table.iloc[:, 1:]
    df_count.to_csv(out_prefix+'.count.xls', sep='\t')


def exp_calculate_and_stat(count_lst,  exp_type='tpm', outdir=os.getcwd(), exp_lower=0.1, outlier_limit=10):
    data = list()
    for each in count_lst:
        table = pd.read_table(each, comment="#", header=0, index_col=[0, 5])
        tmp = table.iloc[:, -1]
        tmp.name = os.path.basename(table.columns[-1]).rsplit('.', 2)[1]
        data.append(tmp)
    df = pd.concat(data, axis=1).reset_index()
    total_out = os.path.join(outdir, 'gene_len_count.txt')
    df.to_csv(total_out, header=True, index=None, sep='\t')
    out_prefix = os.path.join(outdir, os.path.basename(count_lst[0]).split('.', 1)[0])
    exp_calculator_with_count(total_out, exp_type=exp_type, out_prefix=out_prefix,
                              mean_read_len=149, upper_count_percent_limit=0.1)
    # filter very low genes and sort by exp
    exp_file = out_prefix + '.{}.xls'.format(exp_type)
    data = pd.read_table(exp_file, header=0, index_col=0)
    data = data[data['100'] >= exp_lower]
    data = data.sort_values(by='100')
    # data = np.log(data+1)
    # save data
    data.to_csv(exp_file[:-3] + 'filtered.sorted.xls', header=True, index=True, sep='\t')
    # plot deviation
    fig, axes = plt.subplots(2, 2, sharex=True)
    describe = data['100'].describe()
    regions = [
        (describe['min'], describe['25%']),
        (describe['25%'], describe['50%']),
        (describe['50%'], describe['75%']),
        (describe['75%'], describe['max']),
    ]
    subplot_titles = (
        'Q1: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[0]),
        'Q2: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[1]),
        'Q3: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[2]),
        'Q4: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[3]),
    )
    errors = data.apply(lambda x: x / data.loc[:, '100'], axis=0)
    errors = ((errors - 1)*100).abs()
    axes = [y for x in axes for y in x]
    for ind, (lower, upper) in enumerate(regions):
        tmp = errors[(data['100'] >= lower) & (data['100'] <= upper)]
        upper_limit = tmp.describe().loc['50%', :].mean() * outlier_limit
        tmp = tmp[tmp.max(axis=1) <= upper_limit]
        sns.boxplot(data=tmp, whis=3, ax=axes[ind],  showfliers=False)
        axes[ind].tick_params(labelsize="small")
        axes[ind].set_title(subplot_titles[ind], fontsize='small')
    plt.savefig('Exp_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()


def exp_saturation_pipeline(bam, annotation, outdir=os.getcwd(), step=5, samtools="samtools", outlier_limit=10,
                            threads=3, pool_size=4, featureCounts="featureCounts", paired=True):
    bam = os.path.abspath(bam)
    down_bam_lst = down_sample_bam(bam, outdir=outdir, step=step, samtools=samtools,
                                   threads=threads, pool_size=pool_size)
    count_lst = feature_count(down_bam_lst, annotation, featureCounts=featureCounts,
                              paired=paired, threads=threads, pool_size=pool_size)
    exp_calculate_and_stat(count_lst, exp_type='tpm', outdir=outdir, exp_lower=0.5, outlier_limit=outlier_limit)


def introduce_command(func):
    import argparse
    import inspect
    import json
    import time
    parser = argparse.ArgumentParser(description=func.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    func_args = inspect.getfullargspec(func)
    arg_names = func_args.args
    arg_defaults = func_args.defaults
    arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
    for arg, value in zip(arg_names, arg_defaults):
        if value == 'None':
            parser.add_argument('-'+arg, required=True, metavar=arg)
        elif type(value) == bool:
            if value:
                parser.add_argument('--'+arg, action="store_false", help='default: True')
            else:
                parser.add_argument('--'+arg, action="store_true", help='default: False')
        elif value is None:
            parser.add_argument('-' + arg, default=value, metavar='Default:' + str(value), )
        else:
            parser.add_argument('-' + arg, default=value, type=type(value), metavar='Default:' + str(value), )
    if func_args.varargs is not None:
        print("warning: *varargs is not supported, and will be neglected! ")
    if func_args.varkw is not None:
        print("warning: **keywords args is not supported, and will be neglected! ")
    args = parser.parse_args().__dict__
    with open("Argument_detail.json", 'w') as f:
        json.dump(args, f, indent=2, sort_keys=True)
    start = time.time()
    func(**args)
    print("total time: {}s".format(time.time() - start))


if __name__ == '__main__':
    introduce_command(exp_saturation_pipeline)

