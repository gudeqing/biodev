import os
from subprocess import check_call
import pandas as pd
from pysam import VariantFile


def filter_vcf(vcf, out, af=0.02, tumour=1):
    cmd = 'bcftools filter '
    cmd += f'-i "FORMAT/AF[{tumour}:0]>={af}" '
    cmd += f'-o {out} '
    cmd += f'{vcf} '
    print(f'filtering vcf by af >= {af}')
    check_call(cmd, shell=True)
    return out


def annovar_annotation(vcf):
    """
    调用annovar注释vcf而已
    :return:
    """
    script_path = os.path.abspath(__file__)
    if os.path.islink(script_path):
        script_path = os.readlink(script_path)
    dirname = os.path.dirname(script_path)
    check_call(f'sh {dirname}/annovar.sh {vcf} > log.txt', shell=True)
    return f'{vcf}.hg19_multianno.vcf'


def process_annovar_txt(infile, comm_trans, hots=None, af=0.02):
    raw = pd.read_csv(infile, header=0, sep=None, engine='python')
    raw['AF'] = raw[raw.columns[-1]].str.split(':', expand=True)[2].astype(float)
    # 提取常用转录本hgvs注释信息
    cts = dict(x.strip().split()[:2] for x in open(comm_trans))
    hgvs_header = ['Transcript', 'NT_Change', 'AA_Change']
    hgvs_list = []
    for each in raw['AAChange_refGene']:
        hgvs_list.append(['', '', ''])
        if len(each) >= 2:
            all_hgvs = each.split(',')
            for hgvs in all_hgvs:
                hgvs_dict = dict(zip(['gene', 'transcript', 'exon', 'chgvs', 'phgvs'], hgvs.split(':')))
                if cts[hgvs_dict['gene']].startswith(hgvs_dict['transcript'].split('.')[0]) or len(all_hgvs) == 1:
                    hgvs_list[-1] = [
                        hgvs_dict['transcript'],
                        hgvs_dict['chgvs'] if 'chgvs' in hgvs_dict else '',
                        hgvs_dict['phgvs'] if 'phgvs' in hgvs_dict else ''
                    ]
                    break
    hgvs_df = pd.DataFrame(hgvs_list, columns=hgvs_header, index=raw.index)
    raw = raw.join(hgvs_df, rsuffix='_new')
    start_cols = 'Chr,Start,End,Ref,Alt,Func_refGene,GeneDetail_refGene,Gene_refGene,ExonicFunc_refGene,Transcript,NT_Change,AA_Change,AF'.split(',')
    order = start_cols + [x for x in raw.columns if x not in start_cols]
    new = raw.loc[:, order]
    new.to_csv(infile[:-3]+'xls', sep='\t', index=False)

    # 过滤
    filtered = new.loc[new['AF'] >= af, :]
    if filtered.shape[0] <= 0:
        print(f'No mutation pass AF >{af} filtering')
    dirname = os.path.dirname(infile)
    basename = 'filtered.' + os.path.basename(infile)[:-3]+'xls'
    out_filtered = os.path.join(dirname, basename)
    filtered.to_csv(out_filtered, sep='\t', index=False)

    # 提取hot信息
    out = os.path.join(dirname, 'detected.hotspot.xlsx')
    if hots:
        hot_df = pd.read_excel(hots, header=0, index_col=0)
        hot_df.set_index('1', inplace=True)

        # 在new信息中添加is_hotspot
        candidates = new['Otherinfo4'] + ':' + new['Otherinfo5'].astype(str) + ':' + new['Otherinfo6'] + ':' + new['Otherinfo7'] + ':' + new['Otherinfo8']
        new['is_hotspot'] = [x in hot_df.index for x in candidates]
        new['OKR_Name'] = [hot_df.loc[x, 'OKR_Name'] if x in hot_df.index else '' for x in candidates]
        start_cols = start_cols + ['OKR_Name', 'is_hotspot']
        order = start_cols + [x for x in new.columns if x not in start_cols]
        new = new[order]
        new.to_csv(infile[:-3] + 'xls', sep='\t', index=False)

        # 提取hotspot
        candidates = filtered['Otherinfo4']+':'+filtered['Otherinfo5'].astype(str)+':'+filtered['Otherinfo6']+':'+filtered['Otherinfo7']+':'+filtered['Otherinfo8']
        filtered['is_hotspot'] = [x in hot_df.index for x in candidates]
        filtered['OKR_Name'] = [hot_df.loc[x, 'OKR_Name'] if x in hot_df.index else '' for x in candidates]
        filtered = filtered[order]
        filtered.to_csv(out_filtered, sep='\t', index=False)
        hits = [x for x in candidates if x in hot_df.index]
        if hits:
            hits_df = hot_df.loc[hits]
            hits_df.to_excel(out)
        else:
            hits_df = None
            print('No hotspot detected!')
    return out


def extract_hots(vcf, hots, id_mode='chr:start:id:ref:alt', out='detected.hotspot'):
    """
    1. 根据chr:start:ref:alt作为id，取vcf和hots的交集，并以hots的格式输出
    :param vcf:
    :param hots:
    :param id_mode:
    :param out: 输出的excel文件名
    :return:
    """
    detected = []
    id_mode_lst = id_mode.split(':')
    with VariantFile(vcf) as f:
        for record in f:
            # if not list(record.filter)[0] == 'PASS':
            #     continue
            site_info = [record.chrom, str(record.pos), record.id, record.ref, record.alts[0]]
            site_dict = dict(zip(['chr', 'start', 'id', 'ref', 'alt'], site_info))
            site_dict['id'] = '.' # 因为目前hotspot中该信息均为"."
            key = ':'.join(site_dict[x] for x in id_mode_lst)
            if key not in detected:
                detected.append(key)
            else:
                print(f'found duplicated mutation:{key}')
    hot_df = pd.read_excel(hots, header=0, index_col=0)
    hot_df.set_index('1', inplace=True)
    hits = [x for x in detected if x in hot_df.index]
    if hits:
        target = hot_df.loc[hits]
        target.to_excel(out+'.xlsx')
    else:
        print('No hotspot detected')
        target = None
    return target


def pipeline(vcf, hots, comm_trans, af=0.02, tumour=1):
    """
    通过爬虫的方式注释进行okr注释和爬取报告
    :return:
    """
    if not os.path.exists(f'{vcf}.hg19_multianno.vcf'):
        annovar_vcf = annovar_annotation(vcf)
    else:
        print('annovar注释结果已存在,不再重新注释')
        annovar_vcf = f'{vcf}.hg19_multianno.vcf'
    dirname = os.path.dirname(annovar_vcf)
    basename = os.path.basename(annovar_vcf)
    out = os.path.join(dirname, 'filtered.'+basename)
    filtered = filter_vcf(annovar_vcf, out, af=af, tumour=tumour)
    detected = process_annovar_txt(annovar_vcf[:-3]+'txt', comm_trans, hots=hots, af=af)
    if os.path.exists(detected):
        detected = pd.read_excel(detected)
        okr_names = detected['OKR_Name']
        mutations = [x for x in okr_names if type(x) == str]
        print(list(mutations))


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pipeline'])

