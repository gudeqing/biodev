import os
from subprocess import check_call
import pandas as pd
from pysam import VariantFile


def filter_vcf(vcf, out, af=0.02, tumour=1):
    cmd = 'bcftools filter '
    cmd += f'-i "FORMAT/AF[{tumour}:0]>={af}" '
    cmd += f'-o {out} '
    cmd += f'{vcf} '
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
    hot_df = pd.read_csv(hots, header=0, index_col=0)
    hot_df.set_index('1', inplace=True)
    hits = [x for x in detected if x in hot_df.index]
    if hits:
        target = hot_df.loc[hits]
        target.to_excel(out+'.xlsx')
    else:
        print('No hotspot detected')
        target = None
    return target


def pipeline(vcf, hots, af=0.02, tumour=1, id_mode='chr:start:id:ref:alt'):
    """
    通过爬虫的方式注释进行okr注释和爬取报告
    :return:
    """
    annovar_vcf = annovar_annotation(vcf)
    dirname = os.path.dirname(annovar_vcf)
    basename = os.path.basename(annovar_vcf)
    out = os.path.join(dirname, 'filtered.'+basename)
    filtered = filter_vcf(annovar_vcf, out, af=af, tumour=tumour)
    out_name = os.path.join(dirname, 'detected.hotspot')
    detected = extract_hots(filtered, hots, id_mode=id_mode, out=out_name)
    okr_names = detected['OKR_Name']
    mutations = [x for x in okr_names if type(x) == str]
    print(list(mutations))


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pipeline'])

