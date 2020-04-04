import os
import re
import json
import gzip
from collections import namedtuple
import pandas as pd
from pysam import VariantFile
from concurrent.futures import ProcessPoolExecutor as Pool


class VCF(object):
    def __init__(self, vcf):
        """
        :param vcf: 本次解析是参考v4.2版的说明进行的
        """
        self.vcf = vcf
        self.open = gzip.open if vcf.endswith('.gz') else open

    def line2dict(self):
        """
        主要目的是把一行信息转换为字典; header信息除外, 不加处理直接返回。
        返回一个生成器，包含两个元素，第一个元素是字典结构，第二个元素是未加处理的行（方便用于过滤）;
        如果当前行是header信息，则第一个元素是第一个‘=’前的字符，也就是名称。

        header = [ 'CHROM', 'POS', 'ID', 'REF', 'ALT',
            'QUAL', 'FILTER', 'INFO', 'FORMAT',
            'sample1', 'sample2', ...
        ]
        """
        header = list()
        with self.open(self.vcf) as f:
            for line in f:
                if line.startswith('##'):
                    yield line.lstrip('#').split('=', 1)[0], line
                    continue
                elif line.startswith('#CHROM'):
                    header = line.strip('#').strip().split('\t')
                    yield None, line
                    continue
                lst = line.strip().split('\t')
                line_dict = dict(zip(header, lst))
                if 'INFO' in line_dict:
                    info_list = line_dict.pop('INFO').split(';')
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
                        values = line_dict.pop(sample).split(':')
                        if sample in line_dict:
                            print(f'name of sample {sample} are conflicted with others! We will rename it by adding _x')
                            sample += '_x'
                        line_dict[sample] = dict(zip(format_lst, values))
                    line_dict['samples'] = header[9:]
                yield line_dict, line

    def get_field_type(self):
        """
        仅仅提取[INFO，FORMAT, ALT]等相关字段定义信息，获取名称和值的类型
        这些信息可用于制作筛选项，后续可以将这些筛选项组合成filter，用于筛选出符合条件的record
        """
        field_def = dict(INFO=dict(), FORMAT=dict(), ALT=dict())
        with VariantFile(self.vcf) as f:
            for record in f.header.records:
                record_type = record.key
                field_name = record.get('ID')
                field_type = record.get('Type')
                if record_type in ['INFO', 'FORMAT', 'ALT']:
                    field_def[record_type][field_name] = field_type
        return field_def

    def stat_range(self):
        """
        要求vcf中的'.'表示无相关注释，而不是其他含义
        rg表示取值范围 fd表示字段名 tp表示字段取值类型
        本次未使用自己写的line2dict获取信息，主要是为了熟悉pysam的VariantFile
        :return:
        """
        field_type = self.get_field_type()
        value_range = dict()
        with VariantFile(self.vcf) as f:
            for record in f:
                for target in ['INFO', 'FORMAT', 'ALT']:
                    data_container = list()
                    store_container = list()
                    if target == 'INFO':
                        data_container.append(record.info)
                        store_container.append(value_range.setdefault('INFO', {}))
                    elif target == 'FORMAT':
                        for sample in record.samples:
                            data_container.append(record.samples[sample])
                            store_container.append(value_range.setdefault('format_'+sample, {}))
                    else:
                        if field_type['ALT']:
                            data_container.append(record.ALT)
                            store_container.append(record.ALT)
                        else:
                            continue

                    for info, info_dict in zip(data_container, store_container):
                        for fd, tp in field_type[target].items():
                            if fd not in info:
                                continue
                            if type(info[fd]) == tuple:
                                values = info[fd]
                            else:
                                values = [info[fd]]
                            if tp == 'Integer' or tp == 'Float':
                                # 一个字段可能对应多个值，如不同基因型对应的不同统计数据
                                # 以位置序号作为key
                                rg_dict = info_dict.setdefault(fd, dict())
                                for ind, x in enumerate(values):
                                    rg = rg_dict.setdefault(ind, [])
                                    try:
                                        # 值可能是'.'
                                        value = int(x) if tp == 'Integer' else float(x)
                                        if not rg:
                                            rg.append(value)
                                            rg.append(value)
                                        else:
                                            # 比较大小，最后获得值的范围
                                            if value < rg[0]:
                                                rg[0] = value
                                            if value > rg[1]:
                                                rg[1] = value
                                    except:
                                        pass
                            elif tp == 'Character' or tp == 'String':
                                rg_dict = info_dict.setdefault(fd, dict())
                                for ind, x in enumerate(values):
                                    rg = rg_dict.setdefault(ind, set())
                                    if type(x) == str:
                                        rg.add(x.replace('\\x3b', ';').replace('\\x3d', '='))
                                    else:
                                        rg.add(x)
                            else:
                                pass
        # Annovar注释后，有些值都是数值，但是header却说是字符串，现进行校正
        for fd, rg_dict in value_range['INFO'].items():
            for k, v in rg_dict.items():
                if type(v) == set:
                    new_values = []
                    digit = True
                    for x in v:
                        if x != '.':
                            try:
                                new_values.append(float(x))
                            except:
                                digit = False
                                break
                    if digit and new_values:
                        value_range['INFO'][fd][k] = [min(new_values), max(new_values)]
        return value_range

    def export_filter_template(self, name='record', e_num_limit=100):
        """
        :param name:
        :param e_num_limit: 如果某个指标的值非连续且有超过100种, 则不对这个指标进行过滤。
        :return:
        """
        stat_info = self.stat_range()
        rule_dict = dict()
        n = 0
        for section, data in stat_info.items():
            if section == 'INFO':
                # pysam 中的记录用的小写
                section = 'info'
            elif section.startswith('format_'):
                sample = section[7:]
                section = 'samples'
            elif section == 'ALT':
                print("目前不针对此项")
                continue
            # python 语句作为过滤规则，到时直接应用即可
            for fd, rg_dict in data.items():
                for ind, rg in rg_dict.items():
                    if len(rg) <= e_num_limit:
                        n += 1
                        rule_id = f'r{n}'
                        if type(rg) == list and len(rg) == 2:
                            if len(rg_dict) > 1:
                                rule_exp = f"{rg[0]} <= float({name}.{section}['{fd}'][{ind}]) <= {rg[1]} "
                            else:
                                rule_exp = f"{rg[0]} <= float({name}.{section}['{fd}']) <= {rg[1]} "
                            rule_dict[rule_id] = rule_exp
                        elif type(rg) == set and len(rg) >= 2:
                            if section == 'samples':
                                if len(rg_dict) > 1:
                                    rule_exp = f"{name}.{section}['{sample}']['{fd}'][{ind}] in {rg}"
                                else:
                                    rule_exp = f"{name}.{section}['{sample}']['{fd}'] in {rg}"
                                rule_dict[rule_id] = rule_exp
                            else:
                                if len(rg_dict) > 1:
                                    rule_exp = f"{name}.{section}['{fd}'][{ind}] in {rg}"
                                else:
                                    rule_exp = f"{name}.{section}['{fd}'] in {rg}"
                                rule_dict[rule_id] = rule_exp
                    else:
                        print(f'Field {fd} has more than {e_num_limit} kinds of category value, and will not be used!')
        with open('filter_rule_candidates', 'w') as f:
            json.dump(rule_dict, f, indent=2)


class DiffVCF():
    def diff_vcf(self, vcf, vcf2, out='diff.vcf', exclude_if="TRANSLATION_IMPACT=synonymous", true_num=0):
        """
        把vcf2解析成字典，然后读vcf时查询该字典, 排除在字典中存在的记录，然后使用info_filter_if信息进一步过滤。
        'CHROM' + 'POS' + 'REF' + 'ALT' + 'GT'作为字典的key
        :param vcf: 需要处理的vcf文件，将从vcf排除出现在vcf2中的突变
        :param vcf2: 被过滤的信息
        :param out: 过滤后的输出vcf文件
        :param exclude_if: 最后的过滤的条件, 目前仅仅支持"=", 例如"TRANSLATION_IMPACT=synonymous;GENE_REGION=Intronic",
            符合N(=true_num)个条件的记录都会被过滤掉.
        :param true_num: 通过过滤条件的次数, 默认为0，表示要求所有条件都要满足.
        :return:
        """
        tmp_dict = dict()
        for line_dict, line in VCF(vcf2).line2dict():
            if line_dict:
                key = line_dict['CHROM']+line_dict['POS']+line_dict['REF']+line_dict['ALT']
                sample = line_dict['samples'][0]
                sample = sample if sample in line_dict else sample+'_x'
                key += line_dict[sample]['GT']
                tmp_dict[key] = line

        with open(out, 'w') as f:
            for line_dict, line in VCF(vcf).line2dict():
                if type(line_dict)==dict:
                    key = line_dict['CHROM'] + line_dict['POS'] + line_dict['REF'] + line_dict['ALT']
                    sample = line_dict['samples'][0]
                    sample = sample if sample in line_dict else sample + '_x'
                    key += line_dict[sample]['GT']
                    if key not in tmp_dict:
                        if exclude_if:
                            filters = exclude_if.split(';')
                            i = 0
                            for filter_exp in filters:
                                k, v = filter_exp.split('=')
                                if k in line_dict['INFO']:
                                    if line_dict['INFO'].get(k) != v:
                                        i += 1
                                elif k in line_dict:
                                    if line_dict[k] != v:
                                        i += 1
                                else:
                                    # 没有找到，也算通过过滤条件
                                    i += 1
                            if true_num == 0:
                                pass_if = i == len(filters)
                            else:
                                pass_if = i >= true_num
                            if pass_if:
                                f.write(line)
                        else:
                            f.write(line)
                else:
                    f.write(line)






