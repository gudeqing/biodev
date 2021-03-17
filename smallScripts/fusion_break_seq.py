import pysam


def reverse_complement(seq):
    """
    :param seq: 输入序列
    :return: 返回reversed的序列
    """
    seq = seq.upper()
    complement = dict(zip(list("ATCG"), list("TAGC")))
    return ''.join(complement[base] if base in complement else base for base in seq[::-1])


def get_fusion_seq(bedpe, genome='/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta', extend=10, out='target.fa'):
    gn = pysam.FastaFile(genome)
    fields = ['c5', 's5', 'e5', 'c3', 's3', 'e3','name', 'score', 's1', 's2']
    with open(bedpe) as f, open(out, 'w') as fw:
        for line in f:
            lst = line.strip().split()
            ld = dict(zip(fields, lst[:10]))
            if int(ld['e5']) - int(ld['s5']) != 1:
                raise Exception('断点坐标怎么不是相差1')
            # 第一个断点时，正链取上游序列
            if ld['s1'] == '+':
                left = gn.fetch(ld['c5'], int(ld['s5']) - extend, int(ld['e5']))
                # print(len(left))
            else:
                left = gn.fetch(ld['c5'], int(ld['s5']), int(ld['e5']) + extend)
                left = reverse_complement(left)

            # 第二个断点时，正链取下游序列
            if ld['s2'] == '+':
                right = gn.fetch(ld['c3'], int(ld['s3']), int(ld['e3']) + extend)
                # print(len(right))
            else:
                right = gn.fetch(ld['c3'], int(ld['s3']) - extend, int(ld['e3']))
                right = reverse_complement(right)

            # print(left, right)
            b1 = ':'.join([lst[0], lst[-2], lst[2]])
            b2 = ':'.join([lst[3], lst[-1], lst[5]])
            fw.write(f'>{ld["name"]} {b1}--{b2}\n')
            fw.write(left+'|'+right+'\n')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), )

