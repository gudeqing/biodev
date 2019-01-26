import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
import shlex
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


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


def exp_calculator_with_count(count_table_file, exp_type='both', out_prefix=None):
    """
    calculate fpkm and tpm based on count table with second column containing gene length.
    :param count_table_file: example:
    -----------
    gene_id gene_length sample1 sample2
    gene1   1001    29  50
    gene2   1300    30  14
    -----------
    :param exp_type: expression type, fpkm, tpm, or 'both'. default:'both'.
    :return: rpkm_dict, tpm_dict
    """
    if not out_prefix:
        out_prefix = count_table_file[:-4]
    if exp_type.lower() not in ['fpkm', 'tpm', 'both']:
        raise Exception('exp_type should be fpkm or tpm or both')
    count_table = pd.read_table(count_table_file, index_col=0, header=0)
    columns = count_table.columns
    gene_len = count_table[columns[0]]
    if gene_len.min() < 11 or gene_len.max() > 200000:
        print("min gene length: ", gene_len.min())
        print("max gene length: ", gene_len.max())
        print('Warning: The minimum gene length or maximum gene length is abnormal!')
    rpkm_dict = dict()
    tpm_dict = dict()
    for sample in columns[1:]:
        # Divide the read counts by the length of each gene in kilobases.
        # This gives you reads per kilobase (RPK)
        rpk = count_table[sample]/gene_len
        # get rpkm/fpkm
        if exp_type == 'fpkm' or exp_type == 'both':
            total_counts = sum(count_table[sample])
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


def exp_saturation(bam, annotation, outdir=os.getcwd(), step=5, samtools="samtools", exp_type='tpm',
                   threads=3, pool_size=4, featureCounts="featureCounts", paired=True):
    bam = os.path.abspath(bam)
    down_bam_lst = down_sample_bam(bam, outdir=outdir, step=step, samtools=samtools,
                                   threads=threads, pool_size=pool_size)
    count_lst = feature_count(down_bam_lst, annotation, featureCounts=featureCounts,
                              paired=paired, threads=threads, pool_size=pool_size)
    data = list()
    for each in count_lst:
        table = pd.read_table(each, comment="#", header=0, index_col=[0, 5])
        tmp = table.iloc[:, -1]
        tmp.name = os.path.basename(table.columns[-1]).rsplit('.', 2)[1] + b'%'
        data.append(tmp)
    df = pd.concat(data, axis=1).reset_index()
    total_out = os.path.join(outdir, 'gene_len_count.txt')
    df.to_csv(total_out, header=True, index=None, sep='\t')
    out_prefix = os.path.join(outdir, 'all')
    exp_calculator_with_count(total_out, exp_type=exp_type, out_prefix=out_prefix)
    # filter very low genes and sort by exp
    exp_file = out_prefix + '.{}.xls'.format(exp_type)
    data = pd.read_table(exp_file, header=0, index_col=0)
    data = data[(data.sum(axis=1) > 0.1) & (data.iloc[:, -1] < 10)]
    data = data.sort_values(by='100%', axis=0)
    data.to_csv(exp_file[:-3] + '.filtered.sorted.xls', header=True, index=True, sep='\t')
    data = data.apply(lambda x: x/df['100%'], axis=0)
    data.boxplot()
    plt.savefig('Exp_distribution.png', dpi=300, bbox_inches='tight')


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
    introduce_command(exp_saturation)

