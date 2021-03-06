import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import pandas as pd
from collections import OrderedDict
import os
from matplotlib import cm
from matplotlib.colors import Normalize
import numpy as np

trinity_ptr = "~/software/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/DifferentialExpression/PtR"


def expression_var(exp_matrix, exp_threshold=0.1, sep='\t', trinity_ptr=trinity_ptr):
    data = pd.read_table(exp_matrix, index_col=0, header=0, sep=sep)
    expressed = data > exp_threshold
    expressed_num = OrderedDict()
    for each in data.columns:
        expressed_num[each] = expressed[each].value_counts()[True]
    expressed_num_series = pd.Series(expressed_num)
    plt.figure(figsize=(11, 8))
    norm = Normalize(expressed_num_series.min(), expressed_num_series.max())
    color_values = [norm(x) for x in expressed_num_series]
    colors = cm.get_cmap('jet')(color_values)
    # print(colors)
    plt.bar(range(expressed_num_series.shape[0]), expressed_num_series, color=colors)
    plt.xticks(range(expressed_num_series.shape[0]), expressed_num_series.index, rotation=90)
    # expressed_num_series.plot(kind='bar', rot=90, color=colors)
    plt.title('Number of genes with expression over {}'.format(exp_threshold))
    plt.savefig('gene_expression_over_{}_bar.pdf'.format(exp_threshold), dpi=300, bbox_inches='tight')
    # trinity
    if not trinity_ptr:
        exit(0)
    cmd = "perl {} ".format(trinity_ptr)
    cmd += "--matrix {} ".format(exp_matrix)
    cmd += " --min_gene_prevalence {} ".format(data.shape[1] / 2)
    cmd += "--min_gene_expr_val {} ".format(exp_threshold)
    cmd += "--sample_cor_matrix "
    cmd += "--sample_cor_scale_limits 0.1,1 "
    cmd += "--boxplot_log2_dist {} ".format(exp_threshold)
    cmd += " --log2 "
    os.system(cmd)


def introduce_command(func):
    import argparse
    import inspect
    import json
    import time
    parser = argparse.ArgumentParser(description=func.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    func_args = inspect.getargspec(func)
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
    if func_args.keywords is not None:
        print("warning: **keywords args is not supported, and will be neglected! ")
    args = parser.parse_args().__dict__
    with open("Argument_detail.json", 'w') as f:
        json.dump(args, f, indent=2, sort_keys=True)
    start = time.time()
    func(**args)
    print("total time: {}s".format(time.time() - start))


if __name__ == '__main__':
    introduce_command(expression_var)

