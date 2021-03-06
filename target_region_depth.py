import pysam
import pandas as pd
from statistics import median_high


def depth(bed, bams:tuple, out, on_base=False, one_based=False):
    """
    skip reads in which any of the following flags are set:
        BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
    for overlap paired reads, count only once!
    如果按照region统计深度，将取region的中位值作为深度代表
    :param bed: bed file, 至少三列， 如果坐标是1-based，请使用--one_based参数
    :param bam: bam files, separate by white spaces
    :param out: output file
    :param on_base: if stat over each base, default to count over each region defined by each line of bed file.
    :return:
    """
    result = dict()
    for bam in bams:
        sample = bam.rsplit('.', 1)[0].rsplit('/', 1)[1].rsplit('_tumor')[0]
        bam = pysam.AlignmentFile(bam)
        with open(bed) as fr:
            for line in fr:
                if line.startswith('track') or line.startswith('#'):
                    continue
                lst = line.strip().split()
                contig = lst[0]
                start = int(lst[1])-1 if one_based else int(lst[1])
                end = int(lst[2])
                cols = bam.pileup(
                    lst[0], start, end,
                    stepper='all',
                    truncate=True,
                    min_base_quality=10,
                    ignore_orphans=False,
                    ignore_overlaps=True,
                )
                depths = []
                for col in cols:
                    # 利用read名去重overlap
                    depth = len(set(col.get_query_names()))
                    row = [contig, col.reference_pos+1, depth]
                    if on_base:
                        result.setdefault(sample, dict())[tuple(row[:2])] = row[2]
                    else:
                        depths.append(depth)
                if not on_base:
                    row = [contig, start+1, end, median_high(depths)]
                    result.setdefault(sample, dict())[tuple(row[:3])] = row[3]

    data = pd.DataFrame(result)
    data.index.names = ['chr', 'start', 'end'] if not on_base else ['chr', 'pos']
    data = data.sort_index()
    if not on_base:
        bed_info = pd.read_table(bed, header=None, comment='#')
        if not one_based:
            bed_info.iloc[:, 1] = bed_info.iloc[:, 1] + 1
        bed_info = bed_info.set_index(list(bed_info.columns[:3]))
        bed_info.index.names = ['chr', 'start', 'end']
        data = bed_info.join(data)

    if out.endswith('xlsx'):
        data.to_excel(out, merge_cells=False)
    else:
        data.to_csv(out, sep='\t')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['depth'])
