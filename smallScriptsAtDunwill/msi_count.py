import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pysam
import pandas as pd
import re
import statistics
import scipy.stats as stats
from collections import Counter


def parse_region(region):
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
                    'chromosome': lst[0],
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
            yield ld


def msi_stat(region, bam_file, out='None', min_reads=20, cutoff=0.2):
    """
    :param region: 每一行表示一个微卫星位点，可以是MSISensorScan的结果
      example1:  chrM    65      1       2       6       439     953     G       CGTCT   TGTGC
      example2:  chr1 78432506 GGCTT 14[A] GCTAG'格式
    :param bam_file: bam文件路径，需要索引
    :param out: 输出文件
    :param min_reads: 支持某个位点的最小reads数量
    说明:
        判定某个MS是否稳定的一个阈值，默认根据deletion和insertion的频率进行过滤，如比例大于1.5，且频率和大于10%
        根据MSIsenor-pro的文献，MSI样本中，deletion为主要特征，
        而根据经验可推测MSS样本应该比较符合对称的分布如正太分，所以使用上述bias判定稳定性比较可靠
    :return: 字典result[repeat_id] = (repeat_num_lst, alt_ratio, ins_ratio, del_ratio)
    """
    bam = pysam.AlignmentFile(bam_file)
    result = dict()
    print('"^":序列比对的起始位置, "|":比对到repeat区域的起始位置, "$"位置表示比对的终止位置, 如果后续还有序列, 则是未比对上的部分')
    for ld in parse_region(region):
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
            # 计算alt时，如果支持的read数超过2%才算有效，比如该位点有500条reads，那么有效的alt需要至少5个reads
            repeat_num_lst = [x for x in repeat_num_lst if repeat_num_lst.count(x) > int(len(repeat_num_lst)*0.01)+1]
            alt_ratio = sum(x != exp_rp_num for x in repeat_num_lst)/len(repeat_num_lst)
            ins_ratio = sum(x > exp_rp_num for x in repeat_num_lst)/len(repeat_num_lst) + 1e-3
            del_ratio = sum(x < exp_rp_num for x in repeat_num_lst)/len(repeat_num_lst)
        else:
            alt_ratio = 0
            ins_ratio = 0
            del_ratio = 0
            print('NO effective reads found for', repeat_id)
        result[repeat_id] = (repeat_num_lst, alt_ratio, ins_ratio, del_ratio)

    if out != 'None':
        def process(stat_result):
            header = [
                'site', 'bias=del_rate/ins_rate',
                'alt_rate', 'total_reads', 'is_unstable',
                'mean_len', 'len_std', 'len_distribution'
            ]
            lines = []
            site_num = len(stat_result)
            unstable_num = 0
            for k, v in stat_result.items():
                read_num = len(v[0])
                if read_num <= min_reads:
                    site_num -= 1
                    continue
                std = statistics.stdev(v[0])
                mean = statistics.mean(v[0])
                del_ins_bias = round(v[3]/v[2], 2)
                distr = dict(Counter(v[0]))
                # expect = int(k.split('[')[1].split(']')[0])
                # 测试发现MSI-H的细胞系中，大部分sites的分布长度还是符合对称分布的, 将其混入其他标准品中时，会导致明显的不对称性
                # 测试wes临床样本发现，blood样本也会出现不对称性，可能是测序深度不够
                threshold = 1.5
                if 0.1 < v[1] <= 0.13:
                    threshold = 2
                elif 0.13 < v[1] <= 0.2:
                    threshold = 1.8
                unstable = 'no'
                if del_ins_bias >= threshold and v[1] >= 0.1:
                    # 根据MSIsenor-pro的文献，MSS样本中，deletion为主要特征，根据经验可推测，正常样本应该比较符合对称的分布如正太分布
                    # 所以根据insertion和deletion的比例判断是否发生了unstable
                    unstable_num += 1
                    unstable = 'yes'
                lines.append([k, round(del_ins_bias, 2), f'{v[1]:.2%}', read_num, unstable, round(mean, 2), round(std,2), distr])
            lines.sort(key=lambda x: (x[1], x[2]), reverse=True)
            lines = [header] + lines
            status = 'MSI:High' if round(unstable_num/site_num, 2) >= cutoff else 'MSI:Stable'
            lines = [['#Summary', f'{unstable_num}|{site_num}', f'{unstable_num/site_num:.2%}', status]] + lines
            return lines

        lines = process(result)
        with open(out, 'w') as f:
            for line in lines:
                f.write('\t'.join(str(x) for x in line)+'\n')

    return result


def paired_msi(region, normal_bam, tumor_bam, out_prefix='result'):
    print('---Normal---')
    n = msi_stat(region, normal_bam)
    print('---Tumor---')
    t = msi_stat(region, tumor_bam)
    # print(n)
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

