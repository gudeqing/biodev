# coding=utf-8
import os
import pandas as pd


class MergeQuant(object):
    def __init__(self, quant_dir=None, method='salmon', t2g=''):
        self.method = method
        if not quant_dir:
            self.quant_dir = os.getcwd()
        else:
            self.quant_dir = quant_dir
        samples = [os.path.basename(x) for x in os.listdir(quant_dir) if os.path.isdir(x)]
        if not samples:
            raise Exception('quant dir is empty!')
        self.samples = samples
        self.t2g = t2g

    def generate_exp_table(self):
        def merge_file(results, target_cols, new_col_names, out):

            # target_cols = ['TPM', 'count']
            # new_col_names = ['tpm', 'count']
    
            for each_col, new_name in zip(target_cols, new_col_names):
                column_list = list()
                for sample, quant in results:
                    print(sample, quant)
                    tmp_col = pd.read_table(quant, index_col=0, header=0)[each_col]
    
                    # tmp_col.name = sample + '_' + new_name
                    tmp_col.name = sample
                    tmp_col.index.name = 'seq_id'
                    column_list.append(tmp_col)
                if not column_list:
                    return
                result_table = pd.concat(column_list, axis=1)
                result_table.to_csv(out + '.' + new_name + '.matrix', sep='\t')
    
        samples = self.samples
    
        join = os.path.join
    
        if self.method.lower() == 'salmon':
    
            isoform_results = [join(self.quant_dir, x, 'quant.sf') for x in samples]
    
            gene_results = [join(self.quant_dir, x, 'quant.genes.sf') for x in samples]
            print(gene_results)
    
            merge_file(zip(samples, isoform_results), ['TPM', 'NumReads'], ['tpm', 'count'],
    
                       join(self.quant_dir, 'transcript'))
    
            merge_file(zip(samples, gene_results), ['TPM', 'NumReads'], ['tpm', 'count'],
    
                       join(self.quant_dir, 'gene'))
    
        elif self.method.lower() == 'rsem':
    
            iso_results = [join(self.quant_dir, x, x + '.isoforms.results') for x in samples]
    
            gene_results = [join(self.quant_dir, x, x + '.genes.results') for x in samples]
    
            merge_file(zip(samples, iso_results), ['TPM', 'expected_count'], ['tpm', 'count'],
    
                       join(self.quant_dir, 'transcript'))
    
            merge_file(zip(samples, gene_results), ['TPM', 'expected_count'], ['tpm', 'count'],
    
                       join(self.quant_dir, 'gene'))
    
            # get mapping rate
    
            cnt_files = [join(self.quant_dir, x, x + '.stat', x + '.cnt') for x in samples]
    
            with open(join(self.quant_dir, 'alignment_rate.txt'), 'w') as f:
    
                f.write('sample\ttotal_reads\taligned_reads\taligned_rate\n')
    
                for s, each in zip(samples, cnt_files):
                    with open(each) as f2:
                        # ['un-alignable', 'alignable', 'too_many_align', 'total']
    
                        tmp_list = f2.readline().strip('\n').split()
    
                    map_num = int(tmp_list[1]) + int(tmp_list[2])
    
                    map_rate = float(map_num) / int(tmp_list[3])
    
                    f.write('{}\t{}\t{}\t{}\n'.format(s, tmp_list[3], map_num, map_rate))
    
        elif self.method.lower() == 'kallisto':
    
            isoform_results = [join(self.quant_dir, x, 'abundance.tsv') for x in samples]
    
            merge_file(zip(samples, isoform_results), ['tpm', 'est_counts'], ['tpm', 'count'],
    
                       join(self.quant_dir, 'transcript'))
    
            # kallisto do not generate gene expression table, thus t2g file will be needed.
            if not self.t2g or not os.path.exists(self.t2g):
                raise Exception('please provide mapping file: transcript mapping to gene')
            t2g_pd = pd.read_table(self.t2g, index_col=0, header=None, usecols=[0, 1], names=['seq_id', 'gene'])
    
            gene_exp_list = list()
    
            gene_count_list = list()
    
            for sample, result in zip(samples, isoform_results):
                iso_table = pd.read_table(result, index_col=0, header=0)
    
                iso_table = pd.concat([iso_table, t2g_pd], axis=1)
    
                tmp_table = iso_table['est_counts'].groupby(iso_table['gene']).sum()
    
                tmp_table.name = sample
    
                tmp_table.index.name = 'seq_id'
    
                gene_count_list.append(tmp_table)
    
                # exp
    
                tmp_table = iso_table['tpm'].groupby(iso_table['gene']).sum()
    
                tmp_table.name = sample
    
                tmp_table.index.name = 'seq_id'
    
                gene_exp_list.append(tmp_table)
    
            gene_counts = pd.concat(gene_count_list, axis=1)
    
            gene_counts.to_csv(os.path.join(self.quant_dir, 'gene.count.matrix'), sep='\t')
    
            gene_exp = pd.concat(gene_exp_list, axis=1)
    
            gene_exp.to_csv(os.path.join(self.quant_dir, 'gene.tpm.matrix'), sep='\t')
    
        else:
    
            pass


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-quant_dir', default=os.getcwd(), help='quant result dir of salmon or rsem or kallisto')
    parser.add_argument('-method', default='salmon', help='salmon or rsem or kallisto')
    parser.add_argument('-t2g', default='', help='kallisto need provide mapping file: transcript mapping to gene')
    args = parser.parse_args()

    MergeQuant(quant_dir=args.quant_dir, method=args.method, t2g=args.t2g).generate_exp_table()



