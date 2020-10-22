import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pysam
import pandas as pd
import re
import statistics
import scipy.stats as stats
from collections import Counter


def count_msi_per_read(region, bam_file, out='None', cutoff=0.3):
    """
    :param region: msisensor2 MSI格式
    :param bam_file:
    :return:
    """
    bam = pysam.AlignmentFile(bam_file)
    rg_lst = []
    with open(region) as f:
        names = [
            'chromosome',
            'location',
            'repeat_unit_length',
            'repeat_unit_binary',
            'repeat_times',
            'left_flank_binary',
            'right_flank_binary',
            'repeat_unit_bases',
            'left_flank_bases',
            'right_flank_bases'
        ]
        for line in f:
            if '[' in line and ']' in line:
                # for this style: chr1 78432506 GGCTT 14[A] GCTAG
                lst = line.strip().split()
                repeat_unit = lst[3].split('[')[1].split(']')[0]
                ld = {
                    'chromosome':lst[0],
                    'location': lst[1],
                    'left_flank_bases': lst[2],
                    'repeat_unit_bases': repeat_unit,
                    'right_flank_bases': lst[4],
                    'repeat_times': lst[3].split('[')[0],
                    'repeat_unit_length': len(repeat_unit)
                }
            else:
                ld = dict(zip(names, line.strip().split()))
                # print(ld)
            rg_lst.append(ld)

    result = dict()
    print('"^":序列比对的起始位置, "|":比对到repeat区域的起始位置, "$"位置表示比对的终止位置, 如果后续还有序列, 则是未比对上的部分')
    for ld in rg_lst:
        start = int(ld['location'])
        repeat = ld['repeat_unit_bases']
        rp_len = int(ld['repeat_unit_length'])
        left = ld['left_flank_bases']
        right = ld['right_flank_bases']
        end = start + rp_len*int(ld['repeat_times'])
        exp_rp_num = int(ld['repeat_times'])
        repeat_id = f'{ld["chromosome"]}:{start}:{left}|${repeat}[{exp_rp_num}]$|{right}'
        repeat_num_lst = []
        read_set = set()
        for r in bam.fetch(ld['chromosome'], start, end):
            # aln_pos = r.get_reference_positions()
            if r.reference_end is None or r.reference_start is None or r.is_secondary or r.is_duplicate:
                continue
            if r.reference_end >= end and r.reference_start <= start and r.query_name not in read_set:
                read_set.add(r.query_name)  # read name相同的reads只分析一次
                # 判断当前read是否包含MSI序列在参考基因组的区域，然后数repeat
                # 如果当前MSI只有部分包含在read末尾，则没法判断MSI个数, 通过上面的判断可以过滤掉那些reads
                # 这种数的方式也决定了MSI不能太长，否则也没有办法在reads中数到完整的MSI
                aln_seq = r.query_alignment_sequence
                full_seq = r.query_sequence
                rp_start_to_end = aln_seq[start-r.reference_start:]  # repeat比对起始到read末尾的序列
                ini_seq1 = rp_start_to_end
                repeat_num = 0
                # read 以重复区域结尾或开头，那么没办法判断是否扩增，如何识别这样的read？未实现
                # 当MSI区域前存在插入时，导致无法匹配到repeat，需修正
                if aln_seq[:start-r.reference_start].endswith(left) or left.endswith(aln_seq[:start-r.reference_start]):
                    # 确定前面的序列为left flank，则基本可以断定MSI发生插入或替换,所以要进行下面的搜索
                    if not rp_start_to_end.startswith(repeat):
                        rp_start_to_end = rp_start_to_end[rp_start_to_end.find(repeat):]
                else:
                    # 这里可能有点奇怪，但可以校正回那些由于deletion的导致的MSI缩短的
                    # 如果repeat前面的序列和左侧翼不匹配，则通过下面的正则匹配重新找到
                    # match = re.search(f'{left}.*?{right}', full_seq)
                    match = re.search(left+f'({repeat})+?'+right, full_seq)
                    if match:
                        # print('xxx', rp_start_to_end)
                        rp_start_to_end = match.group()[len(left):]
                        # print('xxx', rp_start_to_end)

                seq = rp_start_to_end
                while True:
                    if seq.startswith(repeat):
                        repeat_num += 1
                        seq = seq[rp_len:]
                    else:
                        break
                # repeat 区域往前推，看是否还有扩增
                seq2 = aln_seq[:start-r.reference_start]
                ini_seq2 = seq2
                while True:
                    if seq2.endswith(repeat):
                        repeat_num += 1
                        seq2 = seq2[:-rp_len]
                    else:
                        break

                # 前面没有比对上的区域如果还能匹配repeat，也要加进来
                if not seq2:
                    # 搜索read前面没有比对上的区域
                    seq3 = full_seq[:r.query_alignment_start]
                    while True:
                        if seq3.endswith(repeat):
                            repeat_num += 1
                            seq3 = seq3[:-rp_len]
                        else:
                            break
                if not seq:
                    # 搜索read尾部没有比对上的区域
                    seq4 = full_seq[r.query_alignment_end:]
                    while True:
                        if seq4.startswith(repeat):
                            repeat_num += 1
                            seq4 = seq4[rp_len:]
                        else:
                            break
                #
                repeat_num_lst.append(repeat_num)
                print(f'>find {repeat_num} repeats of {repeat_id} in aligned part of read {r.query_name}:')
                print(
                    # repeat_id,
                    # repeat_num,
                    # r.query_name,
                    r.cigarstring,
                    full_seq[:r.query_alignment_start]
                    +'^'+ini_seq2
                    +'|'+ini_seq1
                    +'$'+full_seq[r.query_alignment_end:]
                )
        if len(repeat_num_lst) > 0:
            # 计算alt时，如果支持的read数超过5%才算有效，比如该位点有500条reads，那么有效的alt需要至少5个reads
            repeat_num_lst = [x for x in repeat_num_lst if repeat_num_lst.count(x) > int(len(repeat_num_lst)*0.02)]
            alt_ratio = sum(x != exp_rp_num for x in repeat_num_lst)/len(repeat_num_lst)
            ins_ratio = sum(x > exp_rp_num for x in repeat_num_lst)/len(repeat_num_lst) + 1e-5
            del_ratio = sum(x < exp_rp_num for x in repeat_num_lst)/len(repeat_num_lst)
        else:
            alt_ratio = 0
            ins_ratio = 0
            del_ratio = 0
            print('NO effective reads found for', repeat_id)
        result[repeat_id] = (repeat_num_lst, alt_ratio, ins_ratio, del_ratio)
    if out != 'None':
        with open(out, 'w') as f:
            site_num = len(result)
            unstable_num = 0
            # f.write('site\trepeat_distribution\tmutation_ratio%\n')
            header = ['site', 'bias=del_rate/ins_rate', 'alt_rate', 'mean_len', 'len_std', 'total_reads', 'len_distribution']
            f.write('\t'.join(header)+'\n')
            for k, v in result.items():
                read_num = len(v[0])
                std = statistics.stdev(v[0])
                mean = statistics.mean(v[0])
                # median = statistics.median(v[0])
                # conf_range = stats.t.interval(0.90, read_num-1, loc=mean, scale=std)
                # 假设正太分布，根据均值和方差算(min, max) that contains alpha percent of the distribution
                # conf_range = stats.norm.interval(0.90, loc=mean, scale=std)
                # conf_range = stats.poisson.interval(0.9, v[1], loc=int(median))
                # conf_range = (round(conf_range[0]), round(conf_range[1]))
                # expect = int(k.split('[')[1].split(']')[0])
                # 把参考重复长度当作期望均值，检验本次采用是否和均值相同
                # t, pvalue = stats.ttest_1samp(v[0], expect)
                # f.write(f'{k}\t{Counter(v[0])}\t{v[1]:.2%}\t{mean:.2f}\t{std:.2f}\t{read_num}\t{v[3]/v[2]:.2f}\t{pvalue}\n')
                f.write(f'{k}\t{v[3]/v[2]:.2f}\t{v[1]:.2%}\t{mean:.2f}\t{std:.2f}\t{read_num}\t{Counter(v[0])}\n')
                # if v[1] > cutoff:
                # if not(conf_range[0]<expect<conf_range[1]):
                # if not(conf_range[0]<=expect<=conf_range[1]) or v[1] > 0.55:
                if v[3] > cutoff*v[2]:
                    # 根据MSIsenor-pro的文献，MSS样本中，长度分布应该是左偏的，根据经验可推测，正常样本应该比较符合对称的分布如正太分布
                    # 所以根据insertion和deletion的比例判断是否发生了unstable
                    unstable_num += 1
            print('Summary:', site_num, unstable_num, "{:.2%}".format(unstable_num/site_num))
            f.write(f'Summary\t{unstable_num}(bias>{cutoff})/{site_num}\t{unstable_num/site_num:.2%}\n')
    return result


def run(region, normal_bam, tumor_bam, out_prefix='result'):
    print('---Normal---:"?"前面的序列表示bam中标记的和参考基因组没有比对上的部分, "|"表示比对的起始位置和终止位置')
    n = count_msi_per_read(region, normal_bam)
    print('---Tumor---:"?"前面的序列表示bam中标记的和参考基因组没有比对上的部分, "|"表示比对的起始位置和终止位置')
    t = count_msi_per_read(region, tumor_bam)
    print(n)
    with open(out_prefix+'.txt', 'w') as f:
        header = [
            'MSI', 'norm_mean', 'tumor_mean', 'diff_pvalue',
            'norm_alt_ratio', 'tumor_alt_ratio',
            'n_effect_depth', 't_effect_depth',
            'norm_detail', 'tumor_detail'
        ]
        f.write('\t'.join(header)+'\n')
        if len(n) <= 5:
            width = 7
            height = 17/5*len(n)
        elif len(n) <= 10:
            width = 8
            height = 19 / 5 * len(n)
        else:
            width = 9
            height = 22 / 5 * len(n)
        fig, axes = plt.subplots(nrows=len(n), ncols=2, figsize=(width, height))
        plt.rcParams['axes.titleweight'] = 'bold'
        # axes = [y for x in axes for y in x]
        for i, key in enumerate(n.keys()):
            n_detail = n[key][0]
            t_detail = t[key][0]
            n_mean = round(statistics.mean(n_detail), 2) if n_detail else 0
            t_mean = round(statistics.mean(t_detail), 2) if t_detail else 0
            s, pvalue = stats.mannwhitneyu(n_detail, t_detail,)
            n_alt_ratio = round(n[key][1], 2)
            t_alt_ratio = round(t[key][1], 2)
            lst = [
                key, n_mean, t_mean, pvalue,
                n_alt_ratio, t_alt_ratio,
                len(n_detail), len(t_detail),
                Counter(n_detail), Counter(t_detail)
            ]
            f.write('\t'.join(str(x) for x in lst)+'\n')

            if n_detail and t_detail:
                # plot density
                s1 = pd.Series(n_detail)
                s2 = pd.Series(t_detail)
                ax = s1.plot.kde(label='Normal', ax=axes[i, 0], alpha=0.7)
                s2.plot.kde(label='Tumor', ax=ax, alpha=0.7)
                ax.tick_params(labelsize='small')
                ax.legend(fontsize='small', loc='upper left')
                ax.set_title(key+f' $pvalue={round(pvalue, 4)}$', fontdict={'fontsize':6}, loc='left')
                # plot histogram
                ax = s1.plot.hist(ax=axes[i, 1], label='Normal', alpha=0.6)
                s2.plot.hist(ax=ax, label='Tumor', alpha=0.6)
                ax.tick_params(labelsize='small')
                ax.yaxis.tick_right()
                ax.legend(fontsize='small', loc='upper left')
                ax.set_title(key+f' $pvalue={round(pvalue, 4)}$', fontdict={'fontsize':6}, loc='left')
            else:
                axes[i, 0].set_title(key+' NoData!', fontdict={'fontsize':6}, loc='left')
                axes[i, 1].set_title(key+' NoData!', fontdict={'fontsize':6}, loc='left')

        plt.savefig(f'{out_prefix}.msi_len_distr.pdf', bbox_inches='tight')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())

