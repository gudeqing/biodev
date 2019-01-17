import os
import argparse
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('-bam', nargs='+')
parser.add_argument('-bed', help='ref bed file')
parser.add_argument('-outdir', default='')
parser.add_argument('-existed_result', default=None, nargs='+')
args = parser.parse_args()


def read_distribution(bam, bed, out):
    """
    read_distribution.py is tool of RseQC
    :param bam:
    :param bed:
    :param out:
    :return:
    """
    cmd = 'read_distribution.py -r {} -i {} > {}'.format(bed, bam, out)
    print(cmd)
    os.system(cmd)
    return out


if args.existed_result:
    out_list = args.existed_result
else:
    out_list = list()
    for each in args.bam:
        out_prefix = os.path.join(args.outdir, os.path.basename(each).split('.', 1)[0])
        out_file = out_prefix + '.read_distribution.txt'
        out_list.append(read_distribution(each, args.bed, out_file))

target_info_list = list()
for each in out_list:
    table = pd.read_table(each, index_col=0, skiprows=4, skipfooter=1, sep='\s+')
    target_info = table.loc[:, 'Tag_count']
    target_info.name = os.path.basename(each).split('.', 1)[0]
    target_info_list.append(target_info)
    target_info.plot(kind='pie', radius=1.1, labeldistance=1)
    plt.savefig('{}.png'.format(each[:-4]), dpi=300)
    plt.close()

all_info = pd.concat(target_info_list, axis=1)
total_out = os.path.join(args.outdir, 'read_distribution.all.txt')
all_info.to_csv(total_out, header=True, index=True, sep='\t')
all_info = all_info/all_info.sum()
plt.figure(figsize=(10, 8))
all_info.transpose().plot(kind='bar', stacked=True, )
plt.xticks(fontsize=7, rotation=30)
img_file = os.path.join(args.outdir, 'read_distribution.png')
plt.savefig(img_file, dpi=300, bbox_inches='tight')
plt.close()

# plot using plotly
import plotly
import plotly.graph_objs as go
df3 = all_info.T
data = [go.Bar(x=df3.index, y=df3[x], name=x) for x in df3.columns]
layout = go.Layout(barmode='stack')
fig = go.Figure(data=data, layout=layout)
html_file = os.path.join(args.outdir, 'read_distribution.html')
plotly.offline.plot(fig, filename=html_file, auto_open=False)
