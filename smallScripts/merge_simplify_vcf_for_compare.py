import re
import pandas as pd


def merge(vcfs:tuple, names:tuple, out_prefix='merged'):
    assert len(vcfs) == len(names)
    vcf_lst = []
    for  vcf, name in zip(vcfs, names):
        d = pd.read_table(vcf, comment='#', index_col=[0,1,2,3,4,5])
        d.columns = [name+'.'+x for x in d.columns[:-1]] + [name]
        vcf_lst.append(d)
    # merge
    e = vcf_lst[0]
    for vcf in vcf_lst[1:]:
        e = e.join(vcf, how='outer')
    e.to_csv(f'{out_prefix}.detail.xls', sep='\t')
    # simplify
    # pattern = re.compile(r'AAChange_refGene=.[^,;]+')
    f1 = lambda x:x.split(':')[0].split('=')[1]+':'+x.split(':')[-1] if type(x)==str else None
    f2 = lambda x:x.split(':')[2]+':'+x.split(':')[3] if type(x)==str else None
    info = e[[x for x in e.columns if x.endswith('INFO')]].applymap(f1)
    af_dp = e[names].applymap(f2)
    filter_ = e[[x for x in e.columns if x.endswith('FILTER')]]
    # result = info.join(af_dp).join(filter_)
    af_dp.columns = info.columns
    result = (info+':'+af_dp).join(filter_)
    result.to_csv(f'{out_prefix}.simplified.xls', sep='\t')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['merge'])
