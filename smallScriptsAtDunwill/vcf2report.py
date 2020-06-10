import os
from subprocess import check_call
import pandas as pd
from pysam import VariantFile, FastaFile

"""
基于分析流程结果，要出具word报告，主要需要走2大步骤：
步骤一，在linux上完成：
    1. 先使用annovar注释vcf，得到格式为txt的注释文件
    2. 基于annovar注释结果txt文件，根据AF过滤；
       依据chr:start:ref:alt从hospot表中提取命中的hotspot
       输出文件 detected.hotspot.xlsx
       输出文件 filtered.TNVariantCaller.mutations.*.vcf.hg19_multianno.xls，其中增加报告需要的一些加工信息。
       输出文件 okr_query.list
    3. 提取msi和tmb信息
    4. 整理CNV信息, 要求基因的每个区间log2cnv >= 1
    5. 注释和整理融合信息, 注释并且提取指定基因的结果
    6. germline的突变结果整理，提取指定基因的突变
    
步骤二，在window上完成
    接下来的处理，需要在window上实现，因为需要进行okr爬虫注释hotspot
    1. 根据hotspot结果，结合msi和tmb信息，进行okr爬虫注释，获取pdf和txt格式报告
    2. 生成word报告，生成报告的流程由张惠开发，我这边提供的输入文件：
        a. okr注释出来的txt格式报告，该文件包含用药等信息
        b. 步骤一生成的过滤后并加工的txt文件filtered.TNVariantCaller.mutations.*.vcf.hg19_multianno.xls，
           该文件包含如下几列信息，可用作报告中的snp/indel的输入信息
           Transcript NT_Change AA_Change AF up_start  down_start  OKR_Name  is_hotspot
        c. msi和tmb信息
        d. cnv信息
        f. 融合信息
        e. 其他，样本信息和实验信息，由其他途径提供
        
为了实现一键化：
    1. 待张惠这边完成word报告生成流程开发后，我这边将改报告流程和okr注释流程串联起来，从而增加流程连贯性。
    2. 通过paromiko包模拟服务器登陆和任务提交，完成annovar注释和word报告生成，okr注释在本地完成
    
"""


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


def process_annovar_txt(infile, comm_trans, genome=None, hots=None, af=0.02, bp_num=25):
    raw = pd.read_csv(infile, header=0, sep=None, engine='python')
    raw['AF'] = raw[raw.columns[-1]].str.split(':', expand=True)[2].astype(float)
    up_streams = []
    down_streams = []
    if genome:
        genome = FastaFile(genome)
        for chr_name, start in zip(raw['Chr'], raw['Start']):
            up_pos = start - bp_num if start >= bp_num else 0
            up = genome.fetch(chr_name, up_pos, start+1)
            down = genome.fetch(chr_name, start-1, start+bp_num)
            up_streams.append(up)
            down_streams.append(down)
        raw['up_start'] = up_streams
        raw['down_start'] = down_streams
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
    if 'up_start' in raw.columns:
        start_cols.append('up_start')
        start_cols.append('down_start')
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
        candidates = new['Otherinfo4'] + ':' + new['Otherinfo5'].astype(str) + \
                     ':' + new['Otherinfo6'] + ':' + new['Otherinfo7'] + ':' + new['Otherinfo8']
        new['is_hotspot'] = [x in hot_df.index for x in candidates]
        new['OKR_Name'] = [hot_df.loc[x, 'OKR_Name'] if x in hot_df.index else '' for x in candidates]
        start_cols = start_cols + ['OKR_Name', 'is_hotspot']
        order = start_cols + [x for x in new.columns if x not in start_cols]
        new = new[order]
        new.to_csv(infile[:-3] + 'xls', sep='\t', index=False)

        # 提取hotspot
        candidates = filtered['Otherinfo4']+':'+filtered['Otherinfo5'].astype(str)+\
                     ':'+filtered['Otherinfo6']+':'+filtered['Otherinfo7']+':'+filtered['Otherinfo8']
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
    这个函数暂时不用了
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


def parse_cnr(cnr, out='cnv.final.txt', single_cutoff=1.0, mean_cutoff=0.0):
    with open(cnr) as f, open(out, 'w') as fw:
        # chromosome, start, end, gene, depth, log2, weight
        header = f.readline().strip().split('\t')
        result = dict()
        for line in f:
            lst = line.strip().split('\t')
            log2cn = float(lst[-1])
            gene_info = lst[3]
            if abs(log2cn) >= single_cutoff and len(gene_info.split(';')) > 1:
                gene_info_dict = dict(x.split('=') for x in gene_info.split(';') if '=' in x)
                if 'ensembl_gn' in gene_info_dict:
                    gene = gene_info_dict['ensembl_gn']
                    result.setdefault(gene, list())
                    result[gene].append(log2cn)
                else:
                    # print('skip NO gene info found line:', line)
                    pass

        fw.write('gene\tmeanLog2CN\n')
        for gene, log2cn_lst in result.items():
            mean_cn = sum(log2cn_lst)/len(log2cn_lst)
            if abs(mean_cn) >= mean_cutoff:
                fw.write(f'{gene}\t{mean_cn}\n')
    return out


def annotate_sv(vcf):
    cmd = f"""bcftools filter -i 'FILTER="PASS" & INFO/PRECISE=1' {vcf} -o filtered.{vcf} """
    cmd2 = f"/nfs2/software/SVAFotate/svafotate.py -i filtered.{vcf} -o svafotate.filtered.{vcf} -ccdg /nfs2/software/SVAFotate/ccdg_sv_afs.bed.gz --gnomad /nfs2/software/SVAFotate/gnomad_sv_afs.bed.gz"
    cmd3 = f"java -Xmx4g -jar /nfs2/software/snpEff/snpEff.jar -v hg19 svafotate.filtered.{vcf} > snpeff.svafotate.filtered.{vcf}"
    cmd4 = f"/nfs2/software/simple_sv_annotation/simple_sv_annotation.py snpeff.svafotate.filtered.{vcf} -o simple.snpeff.svafotate.filtered.{vcf} --gene_list /storage/dqgu/OKR/reportPipe/database/fusion.gene.list --known_fusion_pairs /nfs2/software/simple_sv_annotation/fusion_pairs.txt"
    cmd5 = f"bcftools annotate -x 'INFO/ANN' -o simple.snpeff.svafotate.filtered.{vcf} simple.snpeff.svafotate.filtered.{vcf}"
    cmd6 = f"/nfs2/software/svtools/vcfToBedpe -i simple.snpeff.svafotate.filtered.{vcf} -o simple.snpeff.svafotate.filtered.{vcf}.bedpe"
    descs = ['filtering', 'svafotate pop freq', 'snpeff annotate', 'simplify annotation', 'remove ANN', 'vcf2Bed']
    for each, desc in zip([cmd, cmd2, cmd3, cmd4, cmd5, cmd6], descs):
        print('running:', desc)
        check_call(each, shell=True)
    with open(f"simple.snpeff.svafotate.filtered.{vcf}.bedpe") as f, open('final_fusion.txt', 'w') as fw:
        #0:CHROM_A 1:START_A 2:END_A 3:CHROM_B 4:START_B 5:END_B 6:ID 7:QUAL 8:STRAND_A 9:STRAND_B 10:TYPE 11:FILTER 12:INFO 13:FORMAT 14:T190079D1L18 15:B180518G1L2
        header = f.readline()
        fw.write('fusion_genes\tbreak_positions\ttype\n')
        for line in f:
            lst = line.strip().split('\t')
            if lst[11].lower() not in ['pass']:
                continue
            # SIMPLE_ANN=DEL|EXON_DEL|GNA14|NM_004297.3|Exon1-2del|3,DEL|EXON_DEL|GNA14-AS1|NR_121184.1|Exon2-5del|3;SV_HIGHEST_TIER=3
            if 'simple_ANN=' in lst[12]:
                annot = lst[12].split('simple_ANN=')[1].split(',')
                g1 = annot[0].split('|')[2]
                g2 = annot[1].split('|')[2]
                break_pos = f'{lst[0]}:{(int(lst[1])+int(lst[2]))//2}-{lst[3]}:{(int(lst[4])+int(lst[5]))//2}'
                fw.write(f'{g1}--{g2}\t{break_pos}\t{lst[10]}\n')
            else:
                pass


def pipeline(vcf, hots, comm_trans, genome=None, af=0.02, tumour=1, msi=0, tmb=0):
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
    detected = process_annovar_txt(annovar_vcf[:-3]+'txt', comm_trans, genome=genome, hots=hots, af=af)
    if os.path.exists(detected):
        detected = pd.read_excel(detected)
        okr_names = detected['OKR_Name']
        mutations = [x for x in okr_names if type(x) == str]
        with open('okr_query.list') as f:
            for each in mutations:
                f.write(each+'\n')
            if msi >= 10:
                mutations.append('Microsatellite instability-High')
                f.write('Microsatellite instability-High\n')
            if tmb >= 10:
                f.write('Tumor Mutational Burden\n')
                mutations.append('Tumor Mutational Burden')
        print(mutations)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pipeline', 'parse_cnr', 'annotate_sv'])

