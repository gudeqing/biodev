# -*- coding: utf-8 -*-
import os
import math
import subprocess
import numpy as np
import pandas as pd
import argparse
from pprint import pprint
# from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor as Pool
import glob
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
__author__ = 'gdq'
matplotlib.style.use('ggplot')


def run_script(codes):
    subprocess.call("Rscript {}".format(codes), shell=True)


class PvalueCorrect(object):
    def multtest_correct(self, p_values, method=3):
        """
        1. Bonferroni. ---> bonferroni
        2. Bonferroni Step-down(Holm) ---> Holm
        3. Benjamini and Hochberg False Discovery Rate ---> BH
        4. FDR Benjamini-Yekutieli --->BY
        :param pvalue_list:
        :param method:
        :return: np.array
        """
        pvalue_list = list(p_values)
        n = len(pvalue_list)
        if method == 1:
            fdr = [eachP*n for eachP in pvalue_list]
        elif method == 2:
            sorted_pvalues = sorted(pvalue_list)
            fdr = [eachP*(n - sorted_pvalues.index(eachP)) for eachP in pvalue_list]
        elif method == 3:
            sorted_pvalues = sorted(pvalue_list)
            fdr = [eachP*n/(sorted_pvalues.index(eachP)+1) for eachP in pvalue_list]
        elif method == 4:
            _, fdr = self.fdr_correction(pvalue_list, alpha=0.05, method='negcorr', is_sorted=False)
        fdr = np.array(fdr)
        fdr[fdr > 1] = 1.
        return fdr

    @staticmethod
    def fdr_correction(pvals, alpha=0.05, method='indep', is_sorted=False):
        """pvalue correction for false discovery rate
        This covers Benjamini/Hochberg for independent or positively correlated and
        Benjamini/Yekutieli for general or negatively correlated tests. Both are
        available in the function multipletests, as method=`fdr_bh`, resp. `fdr_by`.
        Parameters
        ----------
        pvals : array_like
            set of p-values of the individual tests.
        alpha : float
            error rate
        method : {'indep', 'negcorr')
        Returns
        -------
        rejected : array, bool
            True if a hypothesis is rejected, False if not
        pvalue-corrected : array
            pvalues adjusted for multiple hypothesis testing to limit FDR
        Notes
        -----
        If there is prior information on the fraction of true hypothesis, then alpha
        should be set to alpha * m/m_0 where m is the number of tests,
        given by the p-values, and m_0 is an estimate of the true hypothesis.
        (see Benjamini, Krieger and Yekuteli)
        """
        def _ecdf(x):
            """
            no frills empirical cdf used in fdrcorrection
            """
            nobs = len(x)
            return np.arange(1, nobs+1)/float(nobs)

        pvals = np.asarray(pvals)
        if not is_sorted:
            pvals_sortind = np.argsort(pvals)
            pvals_sorted = np.take(pvals, pvals_sortind)
        else:
            pvals_sorted = pvals  # alias

        if method in ['i', 'indep', 'p', 'poscorr']:
            ecdffactor = _ecdf(pvals_sorted)
        elif method in ['n', 'negcorr']:
            cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))   #corrected this
            ecdffactor = _ecdf(pvals_sorted) / cm
        else:
            raise ValueError('only indep and negcorr implemented')
        reject = pvals_sorted <= ecdffactor*alpha
        if reject.any():
            rejectmax = max(np.nonzero(reject)[0])
            reject[:rejectmax] = True

        pvals_corrected_raw = pvals_sorted / ecdffactor
        pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
        del pvals_corrected_raw
        pvals_corrected[pvals_corrected > 1] = 1
        if not is_sorted:
            pvals_corrected_ = np.empty_like(pvals_corrected)
            pvals_corrected_[pvals_sortind] = pvals_corrected
            del pvals_corrected
            reject_ = np.empty_like(reject)
            reject_[pvals_sortind] = reject
            return reject_, pvals_corrected_
        else:
            return reject, pvals_corrected


class DiffExpToolbox(PvalueCorrect):
    """
    A toolbox contains several differential analysis tools for RNAseq data.
    Currently, only 3 tools supported: edgeR, DEseq2, DEGseq
    Note: only parts of the functions of each tool are implemented.
    """
    def __init__(self, count_matrix, group_info, cmp_info, exp_matrix=None, exp_type='fpkm',
                 sig_type='pvalue', stat_cutoff=0.05, fc_cutoff=2, padjust_way=3, pool_size=5):
        """
        initial inputs
        :param count_matrix: path of raw count table, '\t' as separator. No duplicated row header !
        :param exp_matrix: path of normalized expression value table,'\t' as separator.
                           if None, the second column of count_matrix must be gene length which
                           will be used to calculate fpkm using edgeR.
        :param group_info: path of group info, file with at least two columns. if no replicate
               exist, just use sample name as group name. Header line starts with '#'.
               --------------------
               #sample group_name  group_name
               s1   group1
               s2   group1  group3
               s3   group2
               s4   group2  group3
               s5   s5
               s6   s6
               --------------------
        :param cmp_info: path of cmp info, file with only two columns. Header line starts with '#'.
               -----------------
               #ctrl    test
               group1   group2
               group2   group3
               s5       s6
               -----------------
        :return: Results will be generated in current directory.
                 tmp/ contain raw results of diff analysis
                 *_vs*.{edgeR, deseq, ...}.diffexp.xls
        """
        self.pool_size = pool_size
        self.count = count_matrix
        self.exp = exp_matrix
        self.exp_type = exp_type
        if sig_type not in ['pvalue', 'padjust']:
            raise Exception('sig_type is not pvalue/padjust')
        self.sig_type = sig_type
        self.fc_cutoff = fc_cutoff
        self.stat_cutoff = stat_cutoff
        self.padjust_way = padjust_way
        self.count_filtered = None
        self.filtered_seqs = []
        # group_info -> dict, group_name as key, list of sample names as values. {group:[s1,s2,]}
        sample_list = list()
        with open(group_info) as f:
            group_dict = dict()
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                tmp_list = line.strip().split()
                sample_list.append(tmp_list[0])
                for g in tmp_list[1:]:
                    group_dict.setdefault(g, list())
                    group_dict[g].append(tmp_list[0])
            for g in group_dict.keys():
                group_dict[g] = sorted(list(set(group_dict[g])))
        self.group_dict = group_dict
        self.samples = sorted(list(set(sample_list)))
        # comparison info -> list. [(ctrl, test), ...]
        with open(cmp_info) as f:
            cmp_list = list()
            error_names = list()
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                tmp_ctrl, tmp_test = line.strip().split()
                if tmp_ctrl not in self.group_dict and tmp_ctrl not in sample_list:
                    error_names.append(tmp_ctrl)
                if tmp_test not in self.group_dict and tmp_test not in sample_list:
                    error_names.append(tmp_test)
                cmp_list.append((tmp_ctrl, tmp_test))
        cmp_list = sorted(list(set(cmp_list)))
        self.cmp_list = cmp_list
        # print sample info
        pprint('group_dict is: ')
        pprint(self.group_dict)
        pprint('comparison list is (ctrl, test): ')
        pprint(self.cmp_list)
        if error_names:
            raise Exception('Each group name of {} is not in {}!'.format(error_names, group_info))

        # check the consistency between group info and count table
        with open(self.count) as f:
            count_samples = f.readline().strip('\n').split('\t')[1:]
        diff = set(self.samples).difference(set(count_samples))
        if diff:
            raise Exception('samples: {} are not contained in count table file'.format(diff))

        # transform count_table and exp_table to python dict
        if self.exp is None:
            with open(self.count) as f:
                count_samples = f.readline().strip('\n').split('\t')[2:]
            if set(self.samples).difference(set(count_samples)):
                raise Exception('gene length column missing')
            self.exp_calculator_with_count(self.count, exp_type=exp_type)
            self.exp = self.count+'.{}.xls'.format(exp_type)
            self.count = str(self.count) + '.count.xls'
        df = pd.read_table(self.count, index_col=0)
        self.count_dicts = df.to_dict('index')
        df = pd.read_table(self.exp, index_col=0)
        self.exp_dicts = df.to_dict('index')
        if sorted(self.count_dicts.keys()) != sorted(self.exp_dicts.keys()):
            cha = set(self.count_dicts.keys()) ^ set(self.exp_dicts.keys())
            raise Exception("The first id column of count table and exp table are different :{} !".format(list(cha)[:10]))

    def filter(self, cutoff=4, passed_number_cutoff=None, output=None, filter_by='count'):
        if output is None:
            output = os.getcwd()
        out_count = os.path.join(output, os.path.basename(self.count) + '_filtered')
        filter_on = self.count if filter_by=='count' else self.exp
        df = pd.read_table(filter_on, index_col=0, header=0)
        sample_num = df.shape[1]
        passed_number_cutoff = int(sample_num / 2)
        ind = df.apply(lambda x: sum(y > cutoff for y in x) >= passed_number_cutoff , axis=1)
        self.filtered_seqs = list(df.index[ind==False])
        df[ind].to_csv(out_count, header=True, index=True, sep='\t')
        self.count_filtered = out_count

    @staticmethod
    def exp_calculator_with_count(count_table_file, exp_type='both'):
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
        if exp_type not in ['fpkm', 'tpm', 'both']:
            raise Exception('exp_type should be fpkm or tpm or both')
        count_table = pd.read_table(count_table_file, index_col=0)
        columns = count_table.columns
        gene_len = count_table[columns[0]]
        if gene_len.min() < 11 or gene_len.max() > 200000:
            print('The minimum gene length and maximum gene length is abnormal!')
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
            df_rpkm.to_csv(count_table_file+'.fpkm.xls', sep='\t')
        if exp_type == 'tpm' or exp_type == 'both':
            df_tpm = pd.DataFrame(tpm_dict)
            df_tpm.to_csv(count_table_file+'.tpm.xls', sep='\t')
        df_count = count_table.iloc[:, 1:]
        df_count.to_csv(count_table_file+'.count.xls', sep='\t')

    def run(self, script_list):
        with Pool(self.pool_size) as pool:
            pool.map(run_script, script_list)

    def __make_result(self, ctrl, test, target_seqs, stat_dict, out_diff_table, out_deg_list):
        if ctrl in self.group_dict:
            ctrl_samples = self.group_dict[ctrl]
        else:
            ctrl_samples = [ctrl]
        if test in self.group_dict:
            test_samples = self.group_dict[test]
        else:
            test_samples = [test]
        with open(out_diff_table, 'w') as f, open(out_deg_list, 'w') as f2:
            count_header = '_count\t'.join(ctrl_samples+test_samples) + '_count'
            tmp_sep = '_' + self.exp_type + '\t'
            exp_header_list = ctrl_samples+test_samples
            if len(ctrl_samples) >= 2:
                exp_header_list.append(ctrl)
            if len(test_samples) >= 2:
                exp_header_list.append(test)
            exp_header = tmp_sep.join(exp_header_list) + '_' + self.exp_type

            f.write('seq_id\t{}\t{}\tlog2fc\tpvalue\tpadjust\tsignificant\tregulate\n'.format(
                count_header, exp_header))
            cmp_samples = ctrl_samples + test_samples
            for seq_id in target_seqs:
                line_list = [seq_id]
                tmp_count_dict = self.count_dicts[seq_id]
                line_list += [tmp_count_dict[x] for x in cmp_samples]
                tmp_exp_dict = self.exp_dicts[seq_id]
                tmp_exp_list = [tmp_exp_dict[x] for x in cmp_samples]
                if len(ctrl_samples) >= 2:
                    tmp_exp_list.append(sum([tmp_exp_dict[x] for x in ctrl_samples])/len(ctrl_samples))
                if len(test_samples) >= 2:
                    tmp_exp_list.append(sum([tmp_exp_dict[x] for x in test_samples])/len(test_samples))
                line_list += tmp_exp_list
                tmp_stat_dict = stat_dict.get(seq_id)
                if tmp_stat_dict:
                    # get fold change
                    tmp_fc = tmp_stat_dict['log2fc']
                    tmp_fc = 0 if tmp_fc != tmp_fc else tmp_fc
                    line_list.append(tmp_fc)
                    # get pvalue
                    tmp_pvalue = tmp_stat_dict['pvalue']
                    tmp_pvalue = 1 if tmp_pvalue != tmp_pvalue else tmp_pvalue
                    line_list.append(tmp_pvalue)
                    # get adjusted pvalue
                    tmp_padjust = tmp_stat_dict['padjust']
                    tmp_padjust = 1 if tmp_padjust != tmp_padjust else tmp_padjust
                    line_list.append(tmp_padjust)
                    # judge significant
                    tmp_stat = tmp_stat_dict[self.sig_type]
                    if abs(tmp_fc) >= math.log(self.fc_cutoff, 2) and tmp_stat <= self.stat_cutoff:
                        line_list.append('yes')
                    else:
                        line_list.append('no')
                    # judge regulate
                    if tmp_fc == 0:
                        reg = 'no change'
                    elif tmp_fc > 0:
                        reg = 'up'
                    else:
                        reg = 'down'
                    line_list.append(reg)
                    # save DEG list
                    if line_list[-2] == 'yes':
                        f2.write(seq_id.rsplit('.', 1)[0] + '\t' + reg + '\n')
                else:
                    line_list += [0, 1, 1, 'untested', 'untested']
                # save
                f.write('\t'.join([str(x) for x in line_list])+'\n')

    def __diff_stat(self, all_stat_dicts, out_stat):
        significant_info = dict()
        total_deg = set()
        cmp_deg_sum = list()
        for each_cmp in all_stat_dicts.keys():
            tmp_table = all_stat_dicts[each_cmp]
            tmp_table = tmp_table[abs(tmp_table['log2fc']) >= math.log(self.fc_cutoff, 2)]
            deg_list = tmp_table[tmp_table[self.sig_type] <= self.stat_cutoff].index.values
            significant_info[each_cmp] = deg_list
            total_deg.update(deg_list)
            cmp_deg_sum.append(len(deg_list))
        with open(out_stat, 'w') as f:
            f.write('seq_id\t{}\tsum\n'.format('\t'.join(all_stat_dicts.keys())))
            sum_info = '\t'.join([str(x) for x in cmp_deg_sum])
            f.write('{}\t{}\t{}\n'.format(len(self.count_dicts.keys()), sum_info, len(total_deg)))
            all_seqs = list(total_deg)+list(set(self.count_dicts.keys())-total_deg)
            for seq_id in all_seqs:
                yes_no_lst = ['yes' if seq_id in significant_info[tmp_cmp] else 'no'
                              for tmp_cmp in all_stat_dicts.keys()]
                f.write('{}\t{}\t{}\n'.format(seq_id, '\t'.join(yes_no_lst), yes_no_lst.count('yes')))

    def DEGseq(self, sep='\t', method='MARS', threshold_kind=5, output=None):
        """
        Differential Analysis with DEGseq. Currently, Only MARS method are Supported.
        :param sep: separator of count_table
        :param method: "LRT", "CTR", "FET", "MARS", "MATR", "FC"
        :param stat_value: pvalue or qvalue cutoff
        :param fold_change: fold change cutoff
        :param threshold_kind: possible kinds are:
            • ‘1’: pValue threshold,
            • ‘2’: zScore threshold,
            • ‘3’: qValue threshold (Benjamini et al. 1995),
            • ‘4’: qValue threshold (Storey et al. 2003),
            • ‘5’: qValue threshold (Storey et al. 2003) and
              Fold-Change threshold on MA-plot are both required (can
              be used only when ‘method="MARS"’).
        :param output: output directory. If None, current directory used.
        :return: Results will be in output directory
        """
        if self.count_filtered is None:
            count_table = self.count
        else:
            count_table = self.count_filtered
        if int(threshold_kind) not in range(6):
            raise NameError("threshold_kind must be one of [0,1,2,3,4,5]")
        with open(count_table) as f:
            header = f.readline().strip('\n').split(sep)

        # using R package
        script_list = list()
        if output is None:
            output = os.getcwd()
        for ctrl, test in self.cmp_list:
            if ctrl in self.group_dict:
                ctrl_ind = ','.join([str(header.index(x)+1) for x in self.group_dict[ctrl]])
            else:
                ctrl_ind = str(header.index(test)+1)
            if test in self.group_dict:
                test_ind = ','.join([str(header.index(x)+1) for x in self.group_dict[test]])
            else:
                test_ind = str(header.index(test)+1)

            script_name = os.path.join(output, 'DEGseq.{}_vs_{}.r'.format(ctrl, test))
            script_list.append(script_name)
            f = open(script_name, 'w')
            f.write('library(DEGseq)\n')
            f.write("## Calculation for {} vs {} \n".format(ctrl, test))
            f.write("ctrl <- readGeneExp(file='{}', geneCol=1, valCol=c({}))\n".format(
                count_table, ctrl_ind))
            f.write('test <- readGeneExp(file="{}", geneCol=1, valCol=c({}))\n'.format(
                count_table, test_ind))
            tmp_output = os.path.join(output, '{}_vs_{}'.format(ctrl, test)+'.degseq.tmp')
            f.write(
                'DEGexp(geneExpMatrix1=test, geneCol1=1, expCol1=c(2:{}), groupLabel1="{}", '
                'geneExpMatrix2=ctrl, geneCol2=1, expCol2=c(2:{}), groupLabel2="{}", '
                'method="{}", rawCount=T, '
                'outputDir="{}")'
                '\n'.format(test_ind.count(',')+2, test,
                            ctrl_ind.count(',')+2, ctrl,
                            method, tmp_output))
            f.close()
        else:
            self.run(script_list)

        # format result
        cmp_result_dirs = [x for x in os.listdir(output) if x.endswith('.degseq.tmp')]
        all_stat_dicts = dict()
        for each in cmp_result_dirs:
            cmp_result = output + '/' + each + '/output_score.txt'
            df = pd.read_table(cmp_result, index_col=0)
            if int(threshold_kind) == 5 or int(threshold_kind) == 4:
                padjust = 'q-value(Storey et al. 2003)'
            elif int(threshold_kind) == 1:
                padjust = 'p-value'
            elif int(threshold_kind) == 3:
                padjust = 'q-value(Benjamini et al. 1995)'
            elif int(threshold_kind) == 2:
                padjust = 'z-score'
            else:
                threshold_kind = 0
            pvalues = df['p-value']
            if threshold_kind:
                stat_df = pd.DataFrame(dict(pvalue=pvalues,
                                            padjust=df[padjust],
                                            log2fc=df['log2(Fold_change) normalized'], ),
                                       index=pvalues.index)
            else:
                correction = self.multtest_correct
                stat_df = pd.DataFrame(dict(pvalue=df['p-value'],
                                            padjust=correction(pvalues, method=self.padjust_way),
                                            log2fc=df['log2(Fold_change) normalized']),
                                       index=pvalues.index)

            target_seqs = list(stat_df[self.sig_type].sort_values().index) + self.filtered_seqs
            ctrl, test = each.split('.degseq.tmp')[0].split('_vs_')
            stat_dict = stat_df.to_dict('index')
            all_stat_dicts['{}_vs_{}'.format(ctrl, test)] = stat_df
            result_table = os.path.join(output, each.split('.tmp')[0] + '.xls')
            result_delist = os.path.join(output, each.split('.tmp')[0] + '.DE.list')
            self.__make_result(ctrl, test, target_seqs, stat_dict, result_table, result_delist)
        else:
            out_stat = os.path.join(output, 'DEGseq_diff_summary.xls')
            self.__diff_stat(all_stat_dicts, out_stat)

    def edgeR(self, dispersion=0.1, output=None, sep='\t'):
        if self.count_filtered is None:
            count_table = self.count
        else:
            count_table = self.count_filtered
        with open(count_table) as f:
            header = f.readline().strip('\n').split(sep)
        if output is None:
            output = os.getcwd()

        script_list = list()
        for ctrl, test in self.cmp_list:
            if ctrl in self.group_dict:
                ctrl_ind = ','.join([str(header.index(x)) for x in self.group_dict[ctrl]])
                ctrl_names = ["'{}'".format(ctrl) for x in self.group_dict[ctrl]]
            else:
                ctrl_ind = str(header.index(test))
                ctrl_names = ["'{}'".format(ctrl)]
            if test in self.group_dict:
                test_ind = ','.join([str(header.index(x)) for x in self.group_dict[test]])
                test_names = ["'{}'".format(test) for x in self.group_dict[test]]
            else:
                test_ind = str(header.index(test))
                test_names = ["'{}'".format(test)]
            script_name = os.path.join(output, 'edgeR.{}_vs_{}.r'.format(ctrl, test))
            script_list.append(script_name)
            f = open(script_name, 'w')
            f.write('library(limma)\n')
            f.write('library(edgeR)\n')
            f.write('counts <- read.table("{}", header=T, row.names=1, sep="{}")\n'.format(
                count_table, sep))
            f.write("## Calculation for {} vs {} \n".format(ctrl, test))
            f.write('tmp_counts <- counts[, c({})]\n'.format(ctrl_ind + ',' + test_ind))
            f.write('tmp_group <- c({})\n'.format(','.join(ctrl_names + test_names)))
            f.write('y <- DGEList(counts=tmp_counts, group=tmp_group)\n')
            f.write('y <- calcNormFactors(y)\n')
            if (',' not in ctrl_ind) or (',' not in test_ind):
                # do sample vs sample if NO replicates
                f.write('result <- exactTest(y, tmp_group, dispersion={})\n'.format(dispersion))
            else:
                f.write('design <- model.matrix(~0+tmp_group)\n')
                f.write('y <- estimateDisp(y, design, robust=F)\n')
                f.write('fit <- glmQLFit(y, design, robust=F)\n')
                f.write('con = makeContrasts(tmp_group{}-tmp_group{}, levels=design)\n'.format(
                    test, ctrl))
                f.write('result <- glmQLFTest(fit, contrast=con)\n')
            f.write('a=topTags(result, n=dim(result$counts)[1], '
                    'adjust.method="BH", sort.by="PValue")\n')
            f.write('write.table(a, "{}/{}_vs_{}.edger.tmp", sep="\\t", row.names=T, '
                    'col.names=NA, quote=FALSE)\n'.format(output, ctrl, test))
            f.close()
        else:
            self.run(script_list)
        # make final report
        cmp_result_dirs = [os.path.join(output, x) for x in os.listdir(output) if x.endswith(
            '.edger.tmp')]
        all_stat_dicts = dict()
        for each in cmp_result_dirs:
            stat_table = pd.read_table(each, index_col=0)
            pvalues = stat_table['PValue']
            df = pd.DataFrame(dict(pvalue=pvalues,
                                   padjust=self.multtest_correct(pvalues, method=self.padjust_way),
                                   log2fc=stat_table['logFC']), index=pvalues.index)
            stat_dict = df.to_dict('index')
            target_seqs = list(df[self.sig_type].sort_values().index) + self.filtered_seqs
            result_table = each.split('.tmp')[0] + '.xls'
            result_delist = each.split('.tmp')[0] + '.DE.list'
            ctrl, test = os.path.basename(each).split('.edger.tmp')[0].split('_vs_')
            all_stat_dicts['{}_vs_{}'.format(ctrl, test)] = df
            self.__make_result(ctrl, test, target_seqs, stat_dict, result_table, result_delist)
        else:
            out_stat = os.path.join(output, 'edgeR_diff_summary.xls')
            self.__diff_stat(all_stat_dicts, out_stat)

    def DESeq2(self, output=None, sep='\t', padjust_way=None):
        """
        :param output: output directory
        :param padjust_way: the method to pvalue correction in R
        :param sep: the separator used in  count table
        Example of using deseq2:
        cs = read.table("egdeR_input.count.xls", header=T, row.names=1, sep="\t")
        colData <- data.frame(row.names=colnames(cs), group=c('c','c','c','e','e','e'))
        dds <- DESeqDataSetFromMatrix(countData=cs, colData=colData, design= ~group)
        dds <- DESeq(dds)
        res <- results(dds, contrast<-c("group", "e", "c"))
        """
        if self.count_filtered is None:
            count_table = self.count
        else:
            count_table = self.count_filtered
        with open(count_table) as f:
            header = f.readline().strip('\n').split(sep)
        if output is None:
            output = os.getcwd()

        script_list = list()
        for ctrl, test in self.cmp_list:
            if ctrl in self.group_dict:
                ctrl_ind = ','.join([str(header.index(x)) for x in self.group_dict[ctrl]])
                ctrl_names = ["'{}'".format(ctrl) for x in self.group_dict[ctrl]]
            else:
                ctrl_ind = str(header.index(test))
                ctrl_names = ["'{}'".format(ctrl)]
            if test in self.group_dict:
                test_ind = ','.join([str(header.index(x)) for x in self.group_dict[test]])
                test_names = ["'{}'".format(test) for x in self.group_dict[test]]
            else:
                test_ind = str(header.index(test))
                test_names = ["'{}'".format(test)]

            script_name = os.path.join(output, 'DESeq2.{}_vs_{}.r'.format(ctrl, test))
            script_list.append(script_name)
            f = open(script_name, 'w')
            f.write('suppressMessages(library(DESeq2))\n')
            f.write('counts <- read.table("{}", header=T, row.names=1, sep="{}")\n'.format(
                count_table, sep))
            f.write("## Calculation for {} vs {} \n".format(ctrl, test))
            f.write('tmp_counts <- counts[, c({})]\n'.format(ctrl_ind + ',' + test_ind))
            f.write('tmp_counts = floor(tmp_counts+0.5)\n')
            f.write('tmp_group <- c({})\n'.format(','.join(ctrl_names + test_names)))
            f.write('colData <- data.frame(row.names=colnames(tmp_counts), group=tmp_group)\n')
            f.write('dds <- DESeqDataSetFromMatrix(countData=tmp_counts, colData=colData, '
                    'design= ~group)\n')
            f.write('dds <- DESeq(dds)\n')
            if padjust_way is None:
                padjust_way = 'BH'
            f.write('res <- results(dds, contrast<-c("group", "{}", "{}"), '
                    'pAdjustMethod="{}")\n'.format(test, ctrl, padjust_way))
            tmp_out_file = os.path.join(output, "{}_vs_{}.deseq2.tmp".format(ctrl, test))
            f.write('write.table(res, file="{}", sep="\\t", quote=F, '
                    'row.names=T, col.names=NA)\n'.format(tmp_out_file))
            f.close()
        else:
            self.run(script_list)

        # make final report
        cmp_result_dirs = [os.path.join(output, x) for x in os.listdir(output) if x.endswith(
            '.deseq2.tmp')]
        all_stat_dicts = dict()
        for each in cmp_result_dirs:
            stat_table = pd.read_table(each, index_col=0)
            pvalues = stat_table['pvalue']
            if padjust_way is None:
                pvalue_correct = self.multtest_correct
                df = pd.DataFrame(dict(pvalue=pvalues,
                                       padjust=pvalue_correct(pvalues, method=self.padjust_way),
                                       log2fc=stat_table['log2FoldChange']), index=pvalues.index)
            else:
                df = pd.DataFrame(dict(pvalue=pvalues,
                                       padjust=stat_table['padj'],
                                       log2fc=stat_table['log2FoldChange']), index=pvalues.index)
            stat_dict = df.to_dict('index')
            target_seqs = list(df[self.sig_type].sort_values().index) + self.filtered_seqs
            result_table = each.split('.tmp')[0] + '.xls'
            result_delist = each.split('.tmp')[0] + '.DE.list'
            ctrl, test = os.path.basename(each).split('.deseq2.tmp')[0].split('_vs_')
            all_stat_dicts['{}_vs_{}'.format(ctrl, test)] = df
            self.__make_result(ctrl, test, target_seqs, stat_dict, result_table, result_delist)
        else:
            out_stat = os.path.join(output, 'DESeq2_diff_summary.xls')
            self.__diff_stat(all_stat_dicts, out_stat)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-count', type=str, required=True,
                        help='path of read count table. Sample name is used as column name. '
                             'Second column must be gene length if "-exp" is None.')
    parser.add_argument('-exp', type=str, default=None,
                        help="path of expression value table, tab as separator."
                             " If None, the second column of count_matrix must be gene length which"
                             " will be used to calculate fpkm or tpm. NOTE: Expression table "
                             "has nothing to do with differential analysis; Only used in report.")
    parser.add_argument('--exp_type', type=str, default="tpm", help='fpkm or tpm. Default: fpkm')
    parser.add_argument('-group', type=str, required=True,
                        help="path of group info file with at least two columns. First column must"
                             " consist of sample names. Other columns consist of group names."
                             "if no replicate exist, just use sample name as group name. "
                             "Header line starts with '#'")
    parser.add_argument('-cmp', type=str, required=True,
                        help="path of comparison info file with only two columns(ctrl vs test)."
                             " Header line starts with '#'")
    parser.add_argument('-method', type=str, default="edgeR", help='DEGseq or edgeR or DESeq2')
    parser.add_argument('--no_filter', default=False, action='store_true',
                        help='Do no filtering. This option will be ignored by default.')
    parser.add_argument('-output', type=str, default=None, help='output directory.')
    parser.add_argument('-pool', type=int, default=5, help='process number for computing')
    parser.add_argument('--plot', default=False, action='store_true', help="do plotting")
    parser.add_argument('-pvalue', type=float, default=0.05, help='p(q)value cutoff. Default: 0.05')
    parser.add_argument('-fc', type=float, default=2.0, help='fold change cutoff. Default: 2.0')
    parser.add_argument('--count_cutoff', type=float, default=4.0,
                        help='count number cutoff for filtering before diff analysis. Default: 4.0')
    parser.add_argument('-filter_by', type=str, default='count', help="filter gene by count or exp")
    parser.add_argument('--passed_number_cutoff', type=int, default=None,
                        help='sample( count > count_cutoff ) number cutoff for filtering before '
                             'diff analysis. Let M=passed_number_cutoff, N=total_sample_number, '
                             'the following event must happen for a gene to be tested: '
                             'Each gene_count in M samples out of N must >= "count_cutoff". '
                             'Default: self-determined')
    parser.add_argument('--degseq_method', type=str, default='MARS',
                        help='method of degseq. Default: MARS')
    parser.add_argument('--degseq_padjust_way', type=int, default=5,
                        help="param of degseq. Integer in [0,1,2,3,4,5]. "
                             " Default: 5 for using qValue (Storey et al. 2003). "
                             "If 0, this option will be ignored; If 2, z-score will be used.")
    parser.add_argument('-sig_type', default="padjust", type=str,
                        help="pvalue or padjust, for diff significance judgement. Default: padjust")
    parser.add_argument('--dispersion', type=float, default=0.1,
                        help='Only used for single sample vs single sample with edgeR.Default: 0.1')
    parser.add_argument('--deseq2_padjust_way', type=str, default=None,
                        help='One of ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",)'
                             'Default: None')
    parser.add_argument('-padjust_way', type=int, default=3,
                        help='Integer in [1,2,3,4]. '
                             'This option is not used if other *_padjust_way specified.'
                             ' 1. Bonferroni. ---> bonferroni;'
                             ' 2. Bonferroni Step-down(Holm) ' '---> Holm;'
                             ' 3. Benjamini and Hochberg False Discovery Rate ---> BH;'
                             ' 4. FDR Benjamini-Yekutieli --->BY'
                             ' Default: 3')

    # ----------------------------------------------------------------------------------------------
    args = parser.parse_args()
    toolbox = DiffExpToolbox(args.count, args.group, args.cmp,
                             exp_matrix=args.exp,
                             exp_type=args.exp_type,
                             sig_type=args.sig_type,
                             fc_cutoff=args.fc,
                             stat_cutoff=args.pvalue,
                             padjust_way=args.padjust_way,
                             pool_size=args.pool)
    if not args.no_filter:
        toolbox.filter(cutoff=args.count_cutoff, output=args.output, filter_by=args.filter_by,
                       passed_number_cutoff=args.passed_number_cutoff)

    if args.method == 'DEGseq':
        toolbox.DEGseq(method=args.degseq_method, threshold_kind=args.degseq_padjust_way,
                       output=args.output, )
    elif args.method == 'edgeR':
        toolbox.edgeR(dispersion=args.dispersion, output=args.output, )
    elif args.method == 'DESeq2':
        toolbox.DESeq2(output=args.output, padjust_way=args.deseq2_padjust_way, )
    else:
        raise Exception('Method {} is not supported'.format(args.method))

    def diff_plot(table):
        df = pd.read_table(table, index_col=0, header=0)
        colors = pd.DataFrame(['gray']*df.shape[0], index=df.index, columns=['color'])
        colors[df['significant'] == 'yes'] = 'red'
        colors[df['regulate'] == "down"] = 'green'
        colors[df['significant'] == 'no'] = 'gray'
        scatter_df = pd.concat([np.log(df.iloc[:, -7]+1), np.log(df.iloc[:, -6]+1), colors], axis=1)
        scatter_df.sort_values(by='color', axis=0, ascending=True, inplace=True)
        x_label, y_label, _ = scatter_df.columns
        scatter_df.plot.scatter(x=x_label, y=y_label, c=scatter_df['color'])
        plt.savefig(os.path.join(args.output, x_label + '_vs_' + y_label + '.scatter.png'), dpi=300)
        plt.close()

        volcano_df = pd.concat([df['log2fc'], -np.log10(df[args.sig_type]), colors], axis=1)
        volcano_df.columns = ['log2fc', '-log10(' + args.sig_type + ')', 'color']
        volcano_df.sort_values(by='color', axis=0, ascending=True, inplace=True)
        x_label2, y_label2, _ = volcano_df.columns
        volcano_df.plot.scatter(x=x_label2, y=y_label2, c=volcano_df['color'])
        plt.savefig(os.path.join(args.output, x_label + '_vs_' + y_label + '.volcano.png'), dpi=300)
        plt.close()

    def density_plot():
        exp = pd.read_table(toolbox.exp, index_col=0, header=0)
        exp = exp[exp.mean(axis=1) >= 0.05]
        # exp_df = np.log(exp+1).dropna()
        exp_df = np.log2(exp+0.01).dropna()
        exp_df.plot(kind="density",)
        plt.xlabel('log2(exp)')
        plt.savefig(os.path.join(args.output, 'exp_based.density.png'), dpi=300)
        plt.close()

        exp = pd.read_table(toolbox.count, index_col=0, header=0)
        exp = exp[exp.mean(axis=1) >= 0.8]
        # exp_df = np.log(exp+1).dropna()
        exp_df = np.log2(exp+0.01).dropna()
        exp_df.plot(kind="density",)
        plt.xlabel('log2(exp)')
        plt.savefig(os.path.join(args.output, 'count_based.density.png'), dpi=300)
        plt.close()

    # plotting
    density_plot()
    if args.plot:
        if args.output is None:
            args.output = os.getcwd()

        results = glob.glob(args.output+'/*_vs_*.{}.xls'.format(args.method.lower()))
        from concurrent.futures import ThreadPoolExecutor as Pool
        with Pool(args.pool) as pool:
            pool.map(diff_plot, results)


