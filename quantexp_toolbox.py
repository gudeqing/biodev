# coding=utf-8
# -*- coding: utf-8 -*-
import os
import subprocess
import shlex
import pandas as pd
import argparse
from multiprocessing import Pool
__author__ = 'gdq'


def run_cmd(cmd):
    # subprocess.check_call should not be used
    # for cannot terminate it when error happens during multiprocessing.
    subprocess.call(cmd)


class QuantExpToolbox(object):
    """
    A toolbox contains several expression quantification tools for RNAseq data.
    salmon: only  quasi-mapping-based mode of Salmon used here.
            http://salmon.readthedocs.io/en/latest/salmon.html
    fastq：a tab separated file with 3 columns, but 3th column is optional. Such as:
    ------------------------------------
    #sample_name left [right]
    sample1 lfq;lfq2;lfq3
    sample2 lfq;lfq2;lfq3   rfq;rfq2;rfq3
    sample3 lfq   rfq
    ------------------------------------

    """
    def __init__(self, transcriptome, fastq, method='salmon', libtype=None,
                 transcript2gene=None, pool=6, threads=10, output=None,
                 salmon=None, rsem=None, kallisto=None,
                 read_len=149, read_len_sd=30,
                 map_tool='bowtie2', map_tool_path=None, samtools=None,
                 ):
        self.method = method
        self.transcriptome = transcriptome
        self.output = os.getcwd() if output is None else output
        self.t2g = transcript2gene
        self.salmon = salmon
        self.rsem = rsem
        self.kallisto = kallisto
        self.pool = pool
        self.threads = threads
        self.fastq = self.parse_fastq(fastq)
        self.libtype = libtype
        self.read_len = read_len
        self.read_len_sd = read_len_sd
        self.map_tool = map_tool
        self.map_tool_path = map_tool_path
        if samtools is None:
            self.samtools = "samtools"
        else:
            self.samtools = samtools

    @staticmethod
    def parse_fastq(fastq):
        fastq_info = dict()
        with open(fastq) as f:
            for line in f:
                if line.startswith('#') or (not line.strip()):
                    pass
                tmp_list = line.strip().split('\t')
                sample, fqs = tmp_list[0], tmp_list[1:]
                fastq_info.setdefault(sample, list())
                read1_list = [x.strip() for x in fqs[0].split(';')]
                fastq_info[sample].append(read1_list)
                if len(fqs) >= 2:
                    read2_list = [x.strip() for x in fqs[1].split(';')]
                    fastq_info[sample].append(read2_list)
        return fastq_info

    def build_index(self, kmer=31):
        if self.method.lower() == 'salmon':
            cmd = '{}/salmon index '.format(self.salmon)
            cmd += '-t {} '.format(self.transcriptome)
            cmd += '-i {}/transcripts_index '.format(self.output)
            cmd += '-k {}'.format(kmer)
            print(cmd)
            subprocess.check_call(shlex.split(cmd))
        elif self.method.lower() == 'kallisto':
            cmd = '{}/kallisto index '.format(self.kallisto)
            cmd += '-i {}/transcripts_index '.format(self.output)
            cmd += '-k {} '.format(kmer)
            cmd += self.transcriptome
            print(cmd)
            subprocess.check_call(shlex.split(cmd))
        elif self.method.lower() == 'rsem':
            if self.map_tool not in ['bowtie', 'bowtie2', 'star']:
                raise Exception('{} is not supported for rsem'.format(self.map_tool))
            cmd = '{}/rsem-prepare-reference '.format(self.rsem)
            if self.t2g is not None:
                t2g = pd.read_table(self.t2g, header=None, index_col=None)
                g2t = t2g.iloc[:, [1, 0]]
                g2t_path = os.path.join(self.output, 'g2t.pair')
                g2t.to_csv(g2t_path, index=False, sep='\t', header=False)
                cmd += "--transcript-to-gene-map {} ".format(g2t_path)
            cmd += "--{} ".format(self.map_tool)
            if self.map_tool_path is not None:
                cmd += "--{}-path {} ".format(self.map_tool, self.map_tool_path)
            cmd += '-p {} '.format(self.threads)
            cmd += '{} '.format(self.transcriptome)
            index_path = '{}/transcripts_index'.format(self.output)
            if not os.path.exists(index_path):
                os.mkdir(index_path)
            cmd += '{}/ref_index'.format(index_path)
            print(cmd)
            subprocess.check_call(shlex.split(cmd))
        else:
            raise Exception('unexpected method: ' + self.method)

    def get_salmon_cmd(self):
        """
        fastq: {sample name: [[fq1,fq2],[fq1,fq2]], }
        :return: cmd list
        """
        cmd_list = list()
        for sample in self.fastq:
            fq_list = self.fastq[sample]
            cmd = '{}/salmon quant '.format(self.salmon)
            cmd += '-i {}/transcripts_index '.format(self.output)
            cmd += '-l A '
            if len(fq_list) == 2:
                cmd += '-1 {} -2 {} '.format(' '.join(fq_list[0]), ' '.join(fq_list[1]))
            else:
                cmd += '-r {} '.format(' '.join(fq_list[0]))
            cmd += '-o {}/{}_quant '.format(self.output, sample)
            cmd += '--gcBias '
            cmd += '-p {} '.format(self.threads)
            if self.t2g is not None:
                cmd += ' -g ' + self.t2g
            cmd_list.append(shlex.split(cmd))
        else:
            return cmd_list

    def get_kallisto_cmd(self):
        cmd_list = list()
        for sample in self.fastq:
            fq_list = self.fastq[sample]
            cmd = '{}/kallisto quant '.format(self.kallisto)
            cmd += '-i {}/transcripts_index '.format(self.output)
            cmd += '-o {}/{}_quant '.format(self.output, sample)
            # cmd += '--plaintext '
            if len(fq_list) == 1:
                cmd += '--single '
                cmd += '-l {} '.format(self.read_len)
                cmd += '-s {} '.format(self.read_len_sd)
            if self.libtype is not None:
                if self.libtype == 'fr':
                    cmd += '--fr-stranded '
                elif self.libtype == 'rf':
                    cmd += '--rf-stranded '
                else:
                    print('the library type {} is invalid and will be ignored'.format(self.libtype))
                    # raise Exception('library type can only be either "fr" or "rf"')
                    pass
            cmd += '-t {} '.format(self.threads)
            if len(fq_list) == 1:
                cmd += ' '.join(fq_list[0])
            else:
                fq_list = list(zip(fq_list[0], fq_list[1]))
                fq_list = [x for y in fq_list for x in y]
                cmd += ' '.join(fq_list)

            cmd_list.append(shlex.split(cmd))
        else:
            return cmd_list

    def get_rsem_cmd(self):
        cmd_list = list()
        for sample in self.fastq:
            fq_list = self.fastq[sample]
            cmd = '{}/rsem-calculate-expression '.format(self.rsem)
            cmd += '-p {} '.format(self.threads)
            if self.libtype is not None:
                if self.libtype == 'fr':
                    cmd += '--forward-prob 1 '
                elif self.libtype == 'rf':
                    cmd += '--forward-prob 0 '
            cmd += '--sort-bam-memory-per-thread 2G '
            cmd += '--estimate-rspd '
            cmd += '--ci-memory 2048 '
            cmd += "--{} ".format(self.map_tool)
            if self.map_tool_path is not None:
                cmd += "--{}-path {} ".format(self.map_tool, self.map_tool_path)
            if len(fq_list) == 1:
                cmd += "--fragment-length-mean {} ".format(self.read_len)
                cmd += "--fragment-length-sd {} ".format(self.read_len_sd)
                cmd += "{} ".format(','.join(fq_list[0]))
            else:
                cmd += '--paired-end '
                cmd += '{} {} '.format(','.join(fq_list[0]), ','.join(fq_list[1]))
            cmd += "{} ".format('{}/transcripts_index/ref_index'.format(self.output))
            out_path = '{}/{}_quant'.format(self.output, sample)
            if not os.path.exists(out_path):
                os.mkdir(out_path)
            cmd += "{}/{}".format(out_path, sample)
            cmd_list.append(shlex.split(cmd))
        else:
            return cmd_list

    def run_quant(self):
        if self.method.lower() == 'salmon':
            cmd_list = self.get_salmon_cmd()
        elif self.method.lower() == 'kallisto':
            cmd_list = self.get_kallisto_cmd()
        elif self.method.lower() == 'rsem':
            cmd_list = self.get_rsem_cmd()
        else:
            raise Exception(self.method + ' is not supported')
        if len(self.fastq) <= self.pool:
            pool = Pool(len(self.fastq))
        else:
            pool = Pool(self.pool)
        pool.map(run_cmd, cmd_list)
        pool.close()
        pool.join()
        self.generate_exp_table()

    def run_sort_bam(self):
        if self.method.lower() == 'rsem':
            samples = sorted(self.fastq.keys())
            bam_list = [os.path.join(self.output, x + '_quant', x + '.transcript.bam') for x in samples]
            cmd_list = list()
            sorted_bam_list = list()
            for each in bam_list:
                out_file = each[:-4] + '.sorted.bam'
                sorted_bam_list.append(out_file)
                cmd = '{} sort --threads {} -o {} {}'.format(self.samtools, self.threads, out_file, each)
                cmd_list.append(shlex.split(cmd))
            bam_list_file = os.path.join(self.output, 'bam.list')
            print(cmd_list)
            with open(bam_list_file, 'w') as f:
                for each in sorted_bam_list:
                    f.write(each + '\n')
            if len(self.fastq) <= self.pool:
                pool = Pool(len(self.fastq))
            else:
                pool = Pool(self.pool)
            pool.map(run_cmd, cmd_list)
            pool.close()
            pool.join()

    def generate_exp_table(self):
        def merge_file(results, target_cols, new_col_names, out):
            # target_cols = ['TPM', 'count']
            # new_col_names = ['tpm', 'count']
            for each_col, new_name in zip(target_cols, new_col_names):
                column_list = list()
                for sample, quant in results:
                    tmp_col = pd.read_table(quant, index_col=0, header=0)[each_col]
                    # tmp_col.name = sample + '_' + new_name
                    tmp_col.name = sample
                    tmp_col.index.name = 'seq_id'
                    column_list.append(tmp_col)
                result_table = pd.concat(column_list, axis=1)
                result_table.to_csv(out+'.'+new_name+'.matrix', sep='\t')

        samples = sorted(self.fastq.keys())
        join = os.path.join
        if self.method.lower() == 'salmon':
            isoform_results = [join(self.output, x+'_quant', 'quant.sf') for x in samples]
            gene_results = [join(self.output, x+'_quant', 'quant.genes.sf') for x in samples]
            merge_file(zip(samples, isoform_results), ['TPM', 'NumReads'], ['tpm', 'count'],
                       join(self.output, 'transcript'))
            merge_file(zip(samples, gene_results), ['TPM', 'NumReads'], ['tpm', 'count'],
                       join(self.output, 'gene'))
        elif self.method.lower() == 'rsem':
            iso_results = [join(self.output, x+'_quant', x+'.isoforms.results') for x in samples]
            gene_results = [join(self.output, x+'_quant', x+'.genes.results') for x in samples]
            merge_file(zip(samples, iso_results), ['TPM', 'expected_count'], ['tpm', 'count'],
                       join(self.output, 'transcript'))
            merge_file(zip(samples, gene_results), ['TPM', 'expected_count'], ['tpm', 'count'],
                       join(self.output, 'gene'))
            # get mapping rate
            cnt_files = [join(self.output, x+'_quant', x+'.stat', x+'.cnt') for x in samples]
            with open(join(self.output, 'alignment_rate.txt'), 'w') as f:
                f.write('sample\ttotal_reads\taligned_reads\taligned_rate\n')
                for s, each in zip(samples, cnt_files):
                    with open(each) as f2:
                        # ['un-alignable', 'alignable', 'too_many_align', 'total']
                        tmp_list = f2.readline().strip('\n').split()
                    map_num = int(tmp_list[1]) + int(tmp_list[2])
                    map_rate = float(map_num)/int(tmp_list[3])
                    f.write('{}\t{}\t{}\t{}\n'.format(s, tmp_list[3], map_num, map_rate))
        elif self.method.lower() == 'kallisto':
            isoform_results = [join(self.output, x+'_quant', 'abundance.tsv') for x in samples]
            merge_file(zip(samples, isoform_results), ['tpm', 'est_counts'], ['tpm', 'count'],
                       join(self.output, 'transcript'))
            # kallisto do not generate gene expression table, thus t2g file will be needed.
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
            gene_counts.to_csv(os.path.join(self.output, 'gene.count.matrix'), sep='\t')
            gene_exp = pd.concat(gene_exp_list, axis=1)
            gene_exp.to_csv(os.path.join(self.output, 'gene.tpm.matrix'), sep='\t')
        else:
            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', type=str, metavar="transcript", required=True,
                        help="transcripts fasta file")
    parser.add_argument('-fq', type=str, metavar="fastq", required=True,
                        help="file contains the fastq path info")
    parser.add_argument('-m', type=str, metavar="method", default='salmon',
                        help="salmon[default], kallisto and rsem are supported.")
    parser.add_argument('-o', type=str, metavar="out_dir", default=None,
                        help="Output directory. Default: current working directory.")
    parser.add_argument('-t2g', type=str, metavar="transcript2gene",
                        help="transcript-->gene, two column: Transcript\tGene. "
                             "But, Salmon also supports gtf or gff input", )
    parser.add_argument('-strand', type=str, metavar="library_type", default=None,
                        help="'fr' or 'rf'. Default: None. Strand-specific protocol type.")
    parser.add_argument('-pool', type=int, metavar="pool_size", default=6,
                        help="Process number for batch analysis of many samples. Default: 6")
    parser.add_argument('-thread', type=str, metavar="thread_number", default=10,
                        help="Threads number for analysis of every single sample.")
    parser.add_argument('-salmon', type=str,
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/Salmon-0.8"
                                ".2_linux_x86_64/bin/",
                        help="where is salmon installed")
    parser.add_argument('-rsem', type=str,
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/RSEM-1.2.31/bin",
                        help="where is rsem installed")
    parser.add_argument('-kallisto', type=str,
                        default="/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/kallisto_linux-v0"
                                ".43.1/",
                        help="where is kallisto installed")
    parser.add_argument('-rl', type=int, metavar="read_length", default=149,
                        help="Only used for single-end mode of kallisto/rsem. "
                             "Estimated average fragment length. Default: 149")
    parser.add_argument('-sd', type=int, metavar="read_length_sd", default=30,
                        help="Only used for single-end mode of kallisto/rsem. "
                             "Use '-sd' along with '-rl. "
                             "Estimated standard deviation of fragment length. Default: 30")
    parser.add_argument('-mapper', type=str, metavar="align_tool", default='bowtie2',
                        help='RSEM argument. Only bowtie, bowtie2 and star are supported. '
                             'Default: bowtie2')
    parser.add_argument('-mapper_path', type=str, metavar="mapper_path",
                        default='/mnt/ilustre/users/sanger-dev/app/bioinfo/align/bowtie2-2.2.9/',
                        help='RSEM argument. Absolute path of alignment tool/mapper')

    args = parser.parse_args()
    toolbox = QuantExpToolbox(transcriptome=args.t,
                              fastq=args.fq,
                              transcript2gene=args.t2g,
                              method=args.m,
                              output=args.o,
                              pool=args.pool,
                              threads=args.thread,
                              libtype=args.strand,
                              salmon=args.salmon,
                              kallisto=args.kallisto,
                              rsem=args.rsem,
                              read_len=args.rl,
                              read_len_sd=args.sd,
                              map_tool=args.mapper,
                              map_tool_path=args.mapper_path,
                              samtools=None,
                              )
    toolbox.build_index()
    toolbox.run_quant()
    toolbox.run_sort_bam()

