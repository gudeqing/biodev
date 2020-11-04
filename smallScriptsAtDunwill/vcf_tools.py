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
        实践表明由于pysam模块可以非常好的处理vcf文件，而且速度快，这里定义的函数可以在vcf不规范时考虑使用
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
                            print(f'name of sample {sample} is conflicted with others! We will rename it by adding _x')
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
        要求vcf中的'.'表示无相关注释，而不是其他含义, 大部分vcf都是遵守这样的规范
        rg表示取值范围 fd表示字段名 tp表示字段取值类型
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
                                # 所有的值都以列表或tuple形式存储
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


def compare_vcf(vcfs:tuple, out, data_fields:tuple=('FORMAT/AF',), often_trans=None, sample_ind=-1,
                comm_ids:tuple=('AAChange_refGeneWithVer',)):
    """
    使用pysam解析vcf，建立结果字典
    mutation唯一id构成：chr:start:ref:alt
    {
        'sample':{'mut1':0.1,'mut2':0.2,'mut3':0.3},
        'sample2':{'mut1':0.12,'mut2':0.22,'mut3':0.32},
        'comm_field': {'mut1':'g:t:c:p','mut2':'g:t:c:p','mut3':'g:t:c:p'},
    }
    :param vcfs: vcf路径信息，空格分开多个vcf
    :param out: 输出文件名，每行以一个突变为单位，第一列为chr:pos:ref:alts构成的突变id，
    :param data_fields: 结果矩阵中的信息对应的字段信息，如果提供多个，将用冒号隔开
    :param often_trans: 常用转录本
    :param sample_ind: 样本索引，即vcf中第几个样本是我们关注的样本，默认提取最后一个样本的信息
    :param comm_ids: 提取其他字段作为每一行的注释，可以是多个字段，提取的信息从输出文件的第二列开始存放。
        默认提取annovar注释中的AAChange_refGene，这种信息是对具体突变的注释，与样本等信息无关，即在每个样本的vcf中都一样。
    :return: 输出out文件，表达矩阵，第一列为chr:pos:ref:alts构成的突变id，第二列开始是突变的描述，剩下列是每个突变对应的量化信息如AF。
    """
    if often_trans:
        often_dict = dict(
            x.strip().split('\t')[:2] for x in open(often_trans) if len(x.strip().split()) > 1
        )
    else:
        often_dict = dict()

    def get_comm_aa_change(changes, often_dict):
        target = None
        if often_dict:
            for each in changes:
                tmp = each.split(':')
                trans = often_dict.get(tmp[0])
                if trans:
                    transcript = tmp[1].split('.')[0]
                    if (trans.split('.')[0] == transcript):
                        target = each
                        break
                else:
                    # print(f'{tmp[0]} is not recorded in common transcripts!')
                    pass
        if target is None:
            target = changes
        return target

    rd = dict()
    for vcf in vcfs:
        with VariantFile(vcf) as f:
            for r in f:
                # mutation id
                mid = ':'.join([str(r.contig), str(r.pos), r.ref, r.alts[0]])
                sample = r.samples.keys()[sample_ind]
                rd.setdefault(sample, dict())
                for each in comm_ids:
                    if each in rd and mid in rd[each]:
                        # 每个id只解析一次
                        continue
                    if each in r.info:
                        rd.setdefault(each, dict())
                        if each == 'AAChange_refGeneWithVer' or each =='AAChange_refGene':
                            cid = get_comm_aa_change(r.info[each], often_dict)
                        else:
                            cid = r.info[each]
                        rd[each][mid] = cid
                values = []
                for each in data_fields:
                    loc, field = each.split('/')
                    if loc == 'FORMAT':
                        if field == 'AF':
                            value = r.samples[sample][field][0]
                        else:
                            value = r.samples[sample][field]
                    else:
                        value = r.info[field]
                    if type(value) == float:
                        value = round(value, 4)
                    values.append(value)
                rd[sample][mid] = '|'.join(str(x) for x in values)

    df = pd.DataFrame(rd)
    df.index.name = 'mutation'
    cids = sorted([x for x in df.columns if x in comm_ids])
    samples = sorted([x for x in df.columns if x not in cids])
    order = cids + samples
    df = df[order]
    df.to_csv(out, sep='\t')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())





