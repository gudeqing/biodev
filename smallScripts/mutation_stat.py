import os
import re
import pandas as pd
from concurrent.futures import ProcessPoolExecutor as Pool
"""
给定一组基因，称之为targets
给定一组样本
计算：
   1. 针对每个样本统计：targets基因的突变情况，相当于grep targets vcf > targets.vcf
   2. 按照1计算所有样本
   3. 得到突变矩阵：行为基因，列为样本
"""


def read_targets(target_file):
    return [x.strip().split()[0] for x in open(target_file)]


def read_sample_metadata(metadata):
    pass


def vcf_generator(vcf):
    header = [
        'CHROM', 'POS', 'ID', 'REF', 'ALT',
        'QUAL', 'FILTER', 'INFO', 'FORMAT',
        'sample1', 'sample2', ...
    ]
    with open(vcf) as f:
        for line in f:
            if line.startswith('##'):
                yield None, line
                continue
            elif line.startswith('#'):
                header = line.strip('#').strip().split('\t')
                yield None, line
                continue
            lst = line.strip().split('\t')
            line_dict = dict(zip(header, lst))
            if 'INFO' in line_dict:
                info_list = line_dict['INFO'].split(';')
                tmp_dict = dict()
                for each in info_list:
                    tmp = [x.strip() for x in each.split('=')]
                    if len(tmp) == 1:
                        tmp_dict[tmp[0]] = 'Flagged'
                    elif len(tmp) == 2:
                        tmp_dict[tmp[0]] = tmp[1]
                    else:
                        raise Exception(f'Parse INFO filed Failed: {each}')
                line_dict['INFO'] = tmp_dict
            if 'FORMAT' in line_dict:
                format_lst = line_dict.pop('FORMAT').split(':')
                for sample in header[9:]:
                    values = line_dict[sample].split(':')
                    if sample in line_dict:
                        print(f'name of sample {sample} are conflicted with other colnames! We will rename it')
                        sample += '_x'
                    line_dict[sample] = dict(zip(format_lst, values))
                line_dict['samples'] = header[9:]
            yield line_dict, line


def filter_paired_vcf(vcf, out='new.vcf'):
    """
    为dna的配对变异分析结果的过滤设计的,当时对比bcftools filter的结果发现，bcftools 似乎并不精确，或者我们对其参数的理解不够精确
    :param vcf:
    :param out:
    :return:
    """
    filter_set = {'strand_bias', 'germline', 'weak_evidence'}
    # filter_set = set()
    tumour_set = set()
    tumour = 'unknown'
    with open(out, 'w') as f:
        for line_dict, line in vcf_generator(vcf):
            if line_dict is None:
                f.write(line)
                continue
            # 推断谁是对照样本
            tumour = line_dict['samples'][0]
            normal = line_dict['samples'][1]
            tumour_gt_sum = sum(float(x) for x in re.split(r'/|\|', line_dict[tumour]['GT']))
            normal_gt_sum = sum(float(x) for x in re.split(r'/|\|', line_dict[normal]['GT']))
            if tumour_gt_sum + normal_gt_sum == 0:
                raise Exception('Found GT of tumour and normal are both zero!')
            if tumour_gt_sum == 0:
                tumour, normal = normal, tumour
            tumour_set.add(tumour)
            if len(tumour_set) != 1:
                print(f'发现当前推断的肿瘤样本和之前的不同: {tumour_set}! while stating')

            tumour_afs = [float(x) for x in line_dict[tumour]['AF'].split(',')]
            normal_afs = [float(x) for x in line_dict[normal]['AF'].split(',')]
            tumour_ad = [float(x) for x in line_dict[tumour]['AD'].split(',')]
            normal_ad = [float(x) for x in line_dict[normal]['AD'].split(',')]

            # 过滤标准如下
            if max(tumour_afs) < 0.05 or \
                    max(tumour_afs) / (max(normal_afs) + 1e-5) < 5 or \
                    sum(tumour_ad) < 20 or \
                    max(tumour_ad[1:]) < 3 or \
                    (set(line_dict['FILTER'].lower().split(';')) & filter_set):
                pass
            else:
                f.write(line)


def target_mutate_status(vcf, target_genes):
    targets = set(read_targets(target_genes))
    status = dict(zip(list(targets), [0] * len(targets)))
    filter_set = {'strand_bias', 'germline', 'weak_evidence'}
    tumour_set = set()
    tumour = 'unknown'
    for line_dict, _ in vcf_generator(vcf):
        if not line_dict:
            continue
        # print('line:', line_dict)
        if 'Gene.refGene' in line_dict['INFO']:
            gene = line_dict['INFO']['Gene.refGene']
        else:
            gene = line_dict['INFO']['Gene_refGene']
        if 'ExonicFunc.refGene' in line_dict['INFO']:
            field = 'ExonicFunc.refGene'
        else:
            field = 'ExonicFunc_refGene'

        # 推断谁是对照样本
        tumour = line_dict['samples'][0]
        normal = line_dict['samples'][1]
        tumour_gt_sum = sum(float(x) for x in re.split(r'/|\|', line_dict[tumour]['GT']))
        normal_gt_sum = sum(float(x) for x in re.split(r'/|\|', line_dict[normal]['GT']))
        if tumour_gt_sum + normal_gt_sum == 0:
            print(gene)
            raise Exception('Found GT of tumour and normal are both zero!')
        if tumour_gt_sum == 0:
            tumour, normal = normal, tumour
        tumour_set.add(tumour)
        if len(tumour_set) != 1:
            print(f'发现当前推断的肿瘤样本和之前的不同: {tumour_set}! while stating {gene}')

        if gene in targets:
            tumour_afs = [float(x) for x in line_dict[tumour]['AF'].split(',')]
            normal_afs = [float(x) for x in line_dict[normal]['AF'].split(',')]
            tumour_ad = [float(x) for x in line_dict[tumour]['AD'].split(',')]
            normal_ad = [float(x) for x in line_dict[normal]['AD'].split(',')]

            # 过滤标准如下
            if max(tumour_afs) < 0.05 or \
                    max(tumour_afs)/(max(normal_afs)+1e-5) < 5 or \
                    sum(tumour_ad) < 20 or \
                    max(tumour_ad[1:]) < 3 or \
                    line_dict['INFO'][field].lower().startswith('synonymous') or \
                    (set(line_dict['FILTER'].lower().split(';')) & filter_set):
                pass
            else:
                # status[gene] = round(max(tumour_afs), 3)
                status[gene] += 1

            # # 以前的过滤标准如下
            # if max(tumour_afs) < 0.05 or \
            #         max(normal_afs) > 0.01 or \
            #         max(tumour_ad[1:]) < 3 or \
            #         sum(tumour_ad) < 30 or \
            #         sum(normal_ad) < 30 or \
            #         line_dict['INFO'][field].lower().startswith('synonymous') or \
            #         (set(line_dict['FILTER'].lower().split(';')) & filter_set):
            #     status[gene] = 0
            # else:
            #     status[gene] = round(max(tumour_afs), 3)

    result = dict()
    result[tumour] = status
    # print(f'Tumour sample is: {tumour}')
    # print(f'Normal sample is: {normal}')
    # print(status)
    return result


def stat_multi_vcf(vcf_lst, target_genes, p_num=5):
    vcfs = list()
    with open(vcf_lst) as f:
        for line in f:
            vcfs.append(line.strip())

    results = list()
    with Pool(p_num) as executor:
        for vcf in vcfs:
            future = executor.submit(target_mutate_status, vcf, target_genes)
            results.append(future)

    merge_result = dict()
    for each in results:
        merge_result.update(each.result())

    df = pd.DataFrame(merge_result)
    df.index.name = 'Genes'
    gene_mutated_ratio = df.apply(lambda x: round(sum(bool(v) for v in x) / df.shape[0], 3), axis=0)
    df.loc['mutated_ratio'] = gene_mutated_ratio
    sample_mutated_ratio = df.apply(lambda x: round(sum(bool(v) for v in x)/df.shape[1], 3), axis=1)
    df['mutated_ratio'] = sample_mutated_ratio
    order = sorted(df.index)
    df = df.loc[order]
    df.to_csv('stats.csv')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['stat_multi_vcf', 'filter_vcf'])

