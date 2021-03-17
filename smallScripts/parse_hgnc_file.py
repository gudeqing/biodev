import pandas as pd


class ParseHGNC(object):
    def __init__(self, data):
        """
        :param data: https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit
        """
        target_cols = [
            'HGNC ID',
            'Approved symbol',
            'Approved name',
            'Previous symbols',
            'Alias symbols',
            'Ensembl gene ID'
        ]
        self.data = pd.read_csv(data, index_col=0, header=0, sep='\t', usecols=target_cols)

    def withdraw_dict(self):
        result = dict()
        for row in self.data.itertuples():
            if '~withdrawn' in row._1 and 'symbol withdrawn' in row._2:
                key = row._1.split('~withdrawn')[0]
                result.setdefault(key, set())
                result[key].add(row._2.split('see')[1].strip())
        return result

    def ensembl2symbols(self):
        """ensembl_id: [approved_name, approved_symbol, alias, alias2, alias3 ...]"""
        result = dict()
        for row in self.data.itertuples():
            if '~withdrawn' in row._1:
                continue
            if pd.isnull(row._5):
                continue
            result.setdefault(row._5, list())
            if not pd.isnull(row._2):
                result[row._5].append(row._2)
            else:
                result[row._5].append('not_found')
            result[row._5].append(row._1)
            if not pd.isnull(row._3):
                result[row._5] += [x.strip() for x in row._3.split(',')]
            if not pd.isnull(row[4]):
                result[row._5] += [x.strip() for x in row[4].split(',')]

        return result

    def symbol2ensembl(self):
        result = dict()
        for k, v in self.ensembl2symbols().items():
            for sym in v[1:]:
                result.setdefault(sym, set())
                result[sym].add(k)
        return result

    def converting(self, query: list or str, out='symbol2ensembl.txt', known_pair=None, symbol2id=True):
        queries = query if type(query) == list else [x.strip().split()[0] for x in open(query)]
        result = list()
        if symbol2id:
            known_dict = dict(x.strip().split()[:2][::-1] for x in open(known_pair) if x.strip()) if known_pair else dict()
        else:
            known_dict = dict(x.strip().split()[:2] for x in open(known_pair) if x.strip()) if known_pair else dict()
        known_dict = {x.lower(): y for x,y in known_dict.items()}
        s2e = self.symbol2ensembl() if symbol2id else self.ensembl2symbols()
        s2e = {x.lower(): y for x,y in s2e.items()}
        with open(out, 'w') as f:
            not_found = []
            for each_ori in queries:
                each = each_ori.lower()
                if each in known_dict:
                    f.write('{}\t{}\n'.format(known_dict[each], each_ori))
                    result.append(known_dict[each])
                else:
                    not_found.append(each_ori)
            if known_pair:
                print('Success to convert {} genes by querying prior known pair'.format(len(result)))

            # find the remained ones
            withdraw_dict = self.withdraw_dict() if symbol2id else dict()
            withdraw_dict = {x.lower(): y for x, y in withdraw_dict.items()}
            failed_ones = []
            for each_ori in not_found:
                each = each_ori.lower()
                if each in s2e:
                    result.append(s2e[each])
                    if len(s2e[each]) > 1:
                        print("{} was found associated with {} genes".format(each_ori, len(s2e[each])))
                    for g in s2e[each]:
                        f.write('{}\t{}\n'.format(g, each_ori))
                elif each in withdraw_dict:
                    print('{} was found in withdraw'.format(each_ori))
                    for new_sym in withdraw_dict[each]:
                        if new_sym in s2e:
                            result.append(s2e[new_sym])
                            for g in s2e[new_sym]:
                                f.write('{}\t{}\n'.format(g, each_ori))
                else:
                    failed_ones.append(each_ori)
                    # print('{} is not found'.format(each))
                    # f.write('{}\t{}\n'.format(each, 'not_found'))
            if not_found:
                print('Success to convert {} genes by query hgnc_custom'.format(len(not_found)-len(failed_ones)))
            if failed_ones:
                print("Failed to query: ")
                print(failed_ones)
        return result


def converting(query, hgnc_custom=None, out='query_result.txt',
               prior_known_pair=None, symbol2id=False):
    """
    converting ensembl id to symbol or reverse
    :param hgnc_custom: https://www.genenames.org/download/custom/， "/nfs2/database/HGNC/custom.txt"
    :param sym: 待查询的列表文件
    :param out: 输出文件名
    :param prior_known_pair: 已经有的ensembl id 和 symbol对应文件, 包含两列; 如提供, 则将优先使用该文件做转换
    :param symbol2id: bool, 如果想把symbol转换为id, 则请用此参数
    """
    hgnc_custom = hgnc_custom if hgnc_custom is not None else 'hgnc.info.txt'
    from urllib.request import urlretrieve
    urlretrieve('https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit', hgnc_custom)
    object = ParseHGNC(hgnc_custom)
    return object.converting(query=query, symbol2id=symbol2id, out=out, known_pair=prior_known_pair)


if __name__ == '__main__':
    from xcmds.xcmds import xcmds
    xcmds(locals(), include=['converting'])








