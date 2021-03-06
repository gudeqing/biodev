import os
import pandas as pd
import sys, time
start = time.time()


class GTF(object):
    def __init__(self, gtf, target_transcripts=None):
        self.gtf = gtf
        self.table = pd.read_csv(self.gtf, sep='\t', comment='#', header=None)
        header = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attrs']
        self.table.columns = header
        # print(self.table.dtypes)
        if target_transcripts:
            # transcripts = set(x.strip().split()[1].split('.')[0] for x in open(comm_trans))
            col9df = pd.DataFrame([self.parse_col9(x) for x in self.table['attrs']], index=self.table.index)
            index = [x.split('.')[0] in target_transcripts if type(x)==str else False for x in col9df['transcript_id']]
            self.table = self.table.loc[index]

    def parse_col9(self, col9: str):
        col9lst = col9.rstrip(';').split(';')
        col9dict = dict(x.strip().replace('"', '').split() for x in col9lst)
        return col9dict

    def single_query(self, chrom, pos, n_nearest=1, tp=None):
        """
        1. 查询的位点是否在某个/某些feature的区域？
        2. 如果不在任何区域，就查找左右最近的feature是谁
        :param chrom: 染色体
        :param pos: 坐标
        :param n_nearest: 指定查找左右最近N个feature，某人只查找一个左右最近的
        :param tp: 限定要查找的feature类型，可以是下面的任何一个，默认所有类型都含括，也即不过滤
                CDS exon gene start_codon stop_codon transcript UTR
        :return: pd.DataFrame, 增加distance列，正数表示左边最近的距离，负数表示右边最近的距离，0表示落在feature区域
        """
        target = self.table[self.table['chr'] == chrom].copy()
        if tp:
            target = target[target['type'] == tp]
        target['pos'] = pos
        locate = target.loc[target['start'] <= pos]
        locate = locate.loc[locate['end'] >= pos]
        if locate.shape[0] == 0:
            # print(f'{chrom}:{pos} locates at no feature region')
            left = target[target['end'] < pos].copy()
            left['distance'] = left['pos'] - left['end']
            left_nearest = left.nsmallest(n_nearest, 'distance', keep='all')
            right = target[target['start'] > pos].copy()
            right['distance'] = right['pos'] - right['start']
            right_nearest = right.nlargest(n_nearest, 'distance', keep='all')
            locate = pd.concat([left_nearest, right_nearest])
        else:
            locate['distance'] = 0
        return locate

    def batch_query(self, chrom_pos_lst, n_nearest=1, tp=None, target_info=None):
        """
        1. 查询的位点是否在某个/某些feature的区域？
        2. 如果不在任何区域，就查找左右最近的feature是谁
        :param chrom_pos_lst: 列表如['chr1:11868', 'chr1:11869']
        :param n_nearest: 指定查找左右最近N个feature，某人只查找一个左右最近的
        :param tp: 限定要查找的feature类型，可以是下面的任何一个，默认所有类型都含括，也即不过滤
                CDS exon gene start_codon stop_codon transcript UTR
        :param target_info: 列表，指定输出感兴趣的attrs，如gene_id, gene_name, transcript_id, exon_number等等
        :return: pd.DataFrame,
            其中distance正数表示查询位点到左边最近的feature距离，负数表示查询位点到右边最近的feature距离，0表示落在feature区域
        """
        results = list()
        for chrom_pos in chrom_pos_lst:
            chrom, pos = chrom_pos.split(':')
            results.append(self.single_query(chrom, int(pos), n_nearest, tp))
        result = pd.concat(results)
        result.set_index(['chr', 'pos', 'distance', 'type'], inplace=True)
        if target_info:
            target_dict = dict()
            col9df = pd.DataFrame([self.parse_col9(x) for x in result['attrs']], index=result.index)
            for col in target_info:
                if col in col9df:
                    target_dict[col] = col9df[col]
                if col in result.columns:
                    target_dict[col] = result[col]
            result = pd.DataFrame(target_dict, index=result.index)
        return result


# result = GTF(sys.argv[1]).batch_query(['chr1:11868', 'chr1:11869', 'chr2:11869'])
result = GTF(sys.argv[1]).batch_query(
    ['chr1:11868', 'chr1:11869', 'chr2:11869',
     'chr10:43114401', 'chr10:43115660', 'chr2:29223566', 'chr6:117325319'],
    tp='exon',
    target_info=['gene_name', 'transcript_id', 'exon_number']
)
print(result)
print('run_time:', round(time.time() - start, 2))




