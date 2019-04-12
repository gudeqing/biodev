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
            'Synonyms',
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

    def query_by_symbol(self, sym: list or str, out='symbol2ensembl.txt'):
        syms = sym if type(sym) == list else [x.strip().split()[0] for x in open(sym)]
        result = list()
        s2e = self.symbol2ensembl()
        withdraw_dict = self.withdraw_dict()
        with open(out, 'w') as f:
            for each in syms:
                if each in s2e:
                    result.append(s2e[each])
                    if len(s2e[each]) > 1:
                        print("{} was found associated with {} genes".format(each, len(s2e[each])))
                    for g in s2e[each]:
                        f.write('{}\t{}\n'.format(each, g))
                elif each in withdraw_dict:
                    print('{} was found in withdraw'.format(each))
                    for new_sym in withdraw_dict[each]:
                        if new_sym in s2e:
                            result.append(s2e[new_sym])
                            for g in s2e[new_sym]:
                                f.write('{}\t{}\n'.format(each, g))
                else:
                    print('{} is not found'.format(each))
                    # f.write('{}\t{}\n'.format(each, 'not_found'))
        return result


if __name__ == '__main__':
    import sys
    obj = ParseHGNC(sys.argv[1])
    obj.query_by_symbol(sys.argv[2], out=sys.argv[3])









