import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib._color_data as mcd
import matplotlib.patches as mpatch
import pandas as pd
import itertools


def get_color_pool(n):
    # https://plot.ly/ipython-notebooks/color-scales/
    return [f'C{x}' for x in range(n)]


def range_lines(bed, header=False, baseline=1, out='lines.pdf', same_height=False):
    """
    bed example:
    chr1	115258671	115258798	NRAS:ENSE00001450282:exon2
    chr1	115256421	115256599	NRAS:ENSE00001450282:exon3
    chr1	115252190	115252349	NRAS:ENSE00001450282:exon4
    chr1	115252173	115252218	chr1:115252173-115252218
    chr1	115252264	115252305	chr1:115252264:115252305
    chr1	115256504	115256546	chr1:115256504:115256546
    chr1	115258729	115258763	chr1:115258729:115258763

    :param bed:
    :param header:
    :param baseline:
    :param out:
    :param same_height:
    :return:
    """
    raw = pd.read_csv(bed, header=0 if header else None, sep=None, engine='python')
    fig, axes = plt.subplots(nrows=baseline)
    if baseline == 1:
        axes = [axes]
    for i, ax in enumerate(axes):
        indexes = [i]+list(range(baseline, raw.shape[0]))
        data = raw.iloc[indexes, :]
        starts = data.iloc[:, 1]
        ends = data.iloc[:, 2]
        names = data.iloc[:, 3]
        colors = get_color_pool(len(starts))
        if not same_height:
            heights = range(len(starts))
        else:
            heights = [0]*len(starts)
        ax.set_yticks(heights)
        ax.set_ylim(-0.1, len(starts)+0.1)
        ax.set_xlim(starts.iloc[0]-1, ends.iloc[0]+1)
        ax.hlines(heights, starts, ends, colors=colors, linewidth=2)
        ax.grid(linestyle='-', linewidth=0.5)
        ax.set_yticklabels([str(x) for x in names], rotation=0, fontsize=6)
        # coord = sorted(pd.concat([starts, ends]))
        # ax.set_xticks(coord)
        # ax.set_xticklabels([str(x) for x in coord], rotation=90, fontsize=6)
    plt.tight_layout()
    plt.savefig(out)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['range_lines'])
