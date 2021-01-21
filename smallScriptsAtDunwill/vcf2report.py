import os
import math
from glob import glob
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


def annot_vcf_by_txt(vcfs:tuple, annot, out, how='inner', mark=None):
    """
    以两个输入文件的前5列信息作为索引，提取hotspot信息
    :param vcfs:
    :param annot:
    :return:
    """
    with open(vcfs[0]) as f:
        for line in f:
            if line.startswith('#CHROM'):
                header = line.split()
                header[-1] = 'FORMAT_Values'
                break

    hotspot = pd.read_table(annot, header=0)
    hotspot.iloc[:, 2] = '.'
    hotspot = hotspot.set_index(list(hotspot.columns[:5]))
    merge = []
    names = []
    for vcf in vcfs:
        name = vcf.split('.')[1]
        names.append(name)
        vcf_df = pd.read_table(vcf, header=None, comment='#')
        vcf_df[2] = '.'  # 抹掉id信息
        # fmt_names = vcf_df[9][0].split(':')
        vcf_df.columns = header
        vcf_df = vcf_df.set_index(list(vcf_df.columns[:5]))
        fmt_df = vcf_df[vcf_df.columns[-1]].str.split(':', expand=True)
        fmt_df = fmt_df.iloc[:, [2, 3]]
        fmt_df.columns = ['AD', 'AF']
        vcf_df = vcf_df.join(fmt_df)

        # print(vcf_df.head())
        # print(hotspot.head())
        result = vcf_df.join(hotspot, how=how, rsuffix='.new')
        result = result.reset_index()
        result.insert(0, 'sample', name)
        merge.append(result)
    all_result = pd.concat(merge)
    all_result.iloc[:, 3] = all_result.iloc[:, 1] + '_' + all_result.iloc[:, 2].astype(str) \
                        + '_' + all_result.iloc[:, 4] + '_' + all_result.iloc[:, 5]
    if mark:
        # mark是一个包含两列的文件，第一列为样本名，第二列为'_'连接的chr_start_ref_alt
        expected = {":".join(x.strip().split()) for x in open(mark)}
        detected = all_result['sample']+':'+ all_result.iloc[:, 3]
        hits = [x in expected for x in detected]
        missed = [x for x in expected if x not in set(detected)]
        if missed:
            print('没有找到', missed)
        else:
            print('所有期望的突变都找到了')
        all_result.insert(1, 'Hit', hits)
    all_result.to_csv(out, sep='\t', index=False)


def filter_vcf_by_af(vcf, af=0.02, tumour=1):
    dirname = os.path.dirname(vcf)
    basename = os.path.basename(vcf)
    out = os.path.join(dirname, 'filtered.' + basename)
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
    if os.path.exists(f'{vcf}.hg19_multianno.vcf'):
        print('annotated result already exists, skip now!')
        return f'{vcf}.hg19_multianno.vcf', f'{vcf}.hg19_multianno.txt'
    annovar = "/nfs2/software/annovar/ANNOVAR_2019.10.24/table_annovar.pl"
    annovardb = "/nfs2/database/annovar/humandb/"
    cmd = f"{annovar} {vcf} {annovardb} -buildver hg19 -out {vcf} -remove "
    g_r_pro = ['refGeneWithVer', 'cytoBand']
    # f_pro = [
    #     'avsnp150', 'cosmic90', 'clinvar_20200316', 'CancerHotspots',
    #     'esp6500siv2_all', '1000g2015aug_all', 'exac03', 'gnomad_exome'
    # ]
    f_pro = [
        'avsnp150', 'snp138NonFlagged',
        'AVENIO', 'EpioneHS',
        'cosmic90', 'Oncomine', 'CGI', 'CBMDB', 'CIVIC', 'DOCM', 'CHASMplus',
        'IntOGen', 'TCGA_PCDM', 'TCGA', 'icgc21', 'CancerHotspots', 'nci60', 'clinvar_20200316',
        'esp6500siv2_all', '1000g2015aug_all', 'exac03', 'gnomad_exome', 'gme', 'cg69', 'hrcr1',
        'intervar_20180118', 'dbnsfp33a'
    ]
    protocols =','.join(g_r_pro+f_pro)
    operations = 'g,r,' + ','.join(['f']*len(f_pro))
    args = ','.join(['-hgvs']*len(g_r_pro+f_pro))
    cmd += f'-protocol {protocols} '
    cmd += f'-operation {operations} '
    cmd += f"--argument '{args}' "
    cmd += "-nastring . -vcfinput --dot2underline --thread 12 --maxgenethread 20 "
    check_call(cmd, shell=True)
    return f'{vcf}.hg19_multianno.vcf', f'{vcf}.hg19_multianno.txt'


def process_annovar_txt(infile, comm_trans, genome=None, hots=None, af=0.02, not_hot_af=0.05, bp_num=25, target_genes=None):
    """
    某些突变有可能不在常用转录本内，而且有多个转录本的注释，这个时候chgvs信息将为空
    :param infile:
    :param comm_trans:
    :param genome:
    :param hots:
    :param af:
    :param not_hot_af:
    :param bp_num:
    :return:
    """
    if type(infile) != str:
        raw = infile
    else:
        raw = pd.read_csv(infile, header=0, sep=None, engine='python')

    if target_genes:
        target_genes = [x.strip() for x in open(target_genes)]
        index = [any(x in target_genes for x in y.split(';')) for y in raw['Gene_refGeneWithVer']]
        raw = raw.loc[index]

    possible_format = raw.iloc[0, -2].split(':')
    af_index = -1
    if 'AF' in possible_format:
        af_index = possible_format.index('AF')
    else:
        possible_format = raw.iloc[0, -3].split(':')
        if 'AF' in possible_format:
            af_index = possible_format.index('AF')
    if af_index > 0:
        raw['AF'] = raw[raw.columns[-1]].str.split(':', expand=True)[af_index].astype(float)
    else:
        try:
            raw['AF'] = raw[raw.columns[-1]].str.split(':', expand=True)[2].astype(float)
        except Exception:
            raw['AF'] = raw[raw.columns[-1]].str.split(':', expand=True)[3].astype(float)

    up_streams = []
    down_streams = []
    if genome:
        genome = FastaFile(genome)
        for chr_name, start in zip(raw['Chr'], raw['Start']):
            up_pos = start - bp_num if start >= bp_num else 0
            up = genome.fetch(chr_name, up_pos, start-1)
            down = genome.fetch(chr_name, start-1, start+bp_num)
            up_streams.append(up)
            down_streams.append(down)
        raw['up_start'] = up_streams
        raw['down_start'] = down_streams
    # 提取常用转录本hgvs注释信息
    cts = dict(x.strip().split()[:2] for x in open(comm_trans))
    hgvs_header = ['Transcript', 'NT_Change', 'AA_Change']
    hgvs_list = []
    aa_changes = raw['AAChange_refGene'] if 'AAChange_refGene' in raw else raw['AAChange_refGeneWithVer']
    for each in aa_changes:
        hgvs_list.append(['', '', ''])
        if len(each) >= 2:
            all_hgvs = each.split(',')
            for hgvs in all_hgvs:
                hgvs_dict = dict(zip(['gene', 'transcript', 'exon', 'chgvs', 'phgvs'], hgvs.split(':')))
                if hgvs_dict['gene'] in cts:
                    if cts[hgvs_dict['gene']].startswith(hgvs_dict['transcript'].split('.')[0]) or len(all_hgvs) == 1:
                        hgvs_list[-1] = [
                            hgvs_dict['transcript'],
                            hgvs_dict['chgvs'] if 'chgvs' in hgvs_dict else '',
                            hgvs_dict['phgvs'] if 'phgvs' in hgvs_dict else ''
                        ]
                        break
                else:
                    print(hgvs_dict['gene'], '没有常用转录本信息')
    hgvs_df = pd.DataFrame(hgvs_list, columns=hgvs_header, index=raw.index)
    raw = raw.join(hgvs_df, rsuffix='_new')
    ver = 'WithVer' if 'Func_refGeneWithVer' in raw else ''
    start_cols = f'Chr,Start,End,Ref,Alt,Func_refGene{ver},GeneDetail_refGene{ver},Gene_refGene{ver},ExonicFunc_refGene{ver},Transcript,NT_Change,AA_Change,AF'.split(',')
    if 'up_start' in raw.columns:
        start_cols.append('up_start')
        start_cols.append('down_start')
    order = start_cols + [x for x in raw.columns if x not in start_cols]
    new = raw.loc[:, order]
    # new是没有过滤的信息，输出加入了新信息的annovar注释结果
    new.to_csv(infile[:-3]+'xls', sep='\t', index=False)
    new.insert(13, 'report', 'yes')

    # 过滤
    filtered = new.loc[new['AF'] >= af, :]
    if filtered.shape[0] <= 0:
        print(f'No mutation pass AF >{af} filtering')
    dirname = os.path.dirname(infile)
    basename = '02.filtered.' + os.path.basename(infile)[:-3]+'xls'
    out_filtered = os.path.join(dirname, basename)
    # 输出过滤后的加入了新信息的annovar注释结果
    filtered.to_csv(out_filtered, sep='\t', index=False)

    # 提取hot信息
    if hots:
        hot_df = pd.read_excel(hots, header=0)
        # print(hot_df.head())
        if 1 in hot_df.columns:
            hot_df.set_index(1, inplace=True)
        elif '1' in hot_df.columns:
            hot_df.set_index('1', inplace=True)
        else:
            hot_df.set_index(hot_df.columns[1], inplace=True)


        # 在new信息中添加is_hotspot
        candidates = new['Otherinfo4'] + ':' + new['Otherinfo5'].astype(str) + \
                     ':' + new['Otherinfo6'] + ':' + new['Otherinfo7'] + ':' + new['Otherinfo8']
        new['is_hotspot'] = [x in hot_df.index for x in candidates]
        new['OKR_Name'] = [hot_df.loc[x, 'OKR_Name'] if x in hot_df.index else '' for x in candidates]
        start_cols = start_cols + ['OKR_Name', 'is_hotspot']
        order = start_cols + [x for x in new.columns if x not in start_cols]
        new = new[order]
        # 在annovar注释的结果中加入OKR_Name和is_hotspot信息后输出
        new.to_csv(infile[:-3] + 'xls', sep='\t', index=False)

        # 按照af过滤结果并输出
        filtered = new.loc[new['AF'] >= af, :]
        # 对于非hotspot采用af>0.05过滤
        inds = filtered['is_hotspot'] | (filtered['AF'] >= not_hot_af)
        filtered = filtered.loc[inds]
        filtered.to_csv(out_filtered, sep='\t', index=False)

        # 提取hotspot，输出命中的hotspot，表格内容以原hotspot为主，即仅仅提取hotspot信息
        # out_name = '01.detected.hotspot.xlsx'
        out_name = '01.hotspot.' + os.path.basename(infile)[:-3] + 'xlsx'
        out = os.path.join(dirname, out_name)
        hits = [x for x in candidates if x in hot_df.index]
        if hits:
            hits_df = hot_df.loc[hits]
            # hits_df['report'] = 'yes'
            hits_df.to_excel(out)
        else:
            hits_df = None
            print('No hotspot detected!')
        return out
    else:
        return out_filtered


def extract_hots(vcf, hots, id_mode='chr:start:id:ref:alt', out='detected.hotspot.xls'):
    """
    这个函数暂时不用了
    1. 根据chr:start:ref:alt作为id，取vcf和hots的交集，并以hots的格式输出
    :param vcf:
    :param hots: excel 文件的hotspot，其中有一列名为1，由chr:start:ref:alt组成，作为唯一id
        如/bid/Hotspot_All/Hotspot_multianno_plus_IVA_OKR_Classification_updateOKR.20200723.xlsx
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
    hot_df = pd.read_excel(hots, header=0)
    hot_df.set_index('1', inplace=True)
    hits = [x for x in detected if x in hot_df.index]
    if hits:
        target = hot_df.loc[hits]
        if out.endswith('xlsx'):
            target.to_excel(out)
        else:
            target.to_csv(out, sep='\t')
    else:
        print('No hotspot detected')
        target = None
    return target


def parse_cnr(cnr, out='cnv.final.txt', amp_cutoff=3, del_cutoff=1.2, okr_targets=None, other_target_genes=None):
    amp_cutoff = math.log2(amp_cutoff / 2)
    del_cutoff = math.log2(del_cutoff / 2)
    okr_cnv_lst = []
    okr_query_lst = []
    if okr_targets:
        okr_cnv_lst = [x.strip() for x in open(okr_targets)]
    if other_target_genes:
        other_target_genes = [x.strip() for x in open(other_target_genes)]
    else:
        other_target_genes = []

    with open(cnr) as f, open(out, 'w') as fw, open(out+'.offtarget', 'w') as fw2:
        # chromosome, start, end, gene, depth, log2, weight
        header = f.readline().strip().split('\t')
        result = dict()
        for line in f:
            lst = line.strip().split('\t')
            log2cn = float(lst[-2])
            gene_info = lst[3]
            if len(gene_info.split(';')) > 1:
                gene_info_dict = dict(x.split('=') for x in gene_info.split(';') if '=' in x)
                if 'ensembl_gn' in gene_info_dict:
                    gene = gene_info_dict['ensembl_gn']
                    result.setdefault(gene, list())
                    result[gene].append(log2cn)
                else:
                    # print('skip NO gene info found line:', line)
                    pass

        fw.write('gene\tmeanLog2CN\tmeanCN\tOKR_Name\treport\n')
        fw2.write('gene\tmeanLog2CN\tmeanCN\tOKR_Name\treport\n')
        for gene, log2cn_lst in result.items():
            mean_cn = sum(log2cn_lst)/len(log2cn_lst)
            exact_cn = round(2**mean_cn*2, 2)
            # if (mean_cn >= amp_cutoff or mean_cn <= del_cutoff):
            if all(x > amp_cutoff for x in log2cn_lst) or all(x <= del_cutoff for x in log2cn_lst):
                # print(log2cn_lst)
                okr_name = f'{gene} amplification' if mean_cn >0 else f'{gene} deletion'
                # if gene == 'ERBB2':
                #     print(okr_name)
                #     print(log2cn_lst)
                if okr_name in okr_cnv_lst:
                    okr_query_lst.append(okr_name)
                    fw.write(f'{gene}\t{mean_cn:.2f}\t{exact_cn}\t{okr_name}\tyes\n')
                elif gene in other_target_genes:
                    fw.write(f'{gene}\t{mean_cn:.2f}\t{exact_cn}\t\tyes\n')
                else:
                    fw2.write(f'{gene}\t{mean_cn:.2f}\t{exact_cn}\t\tno\n')
    print(f'hit {len(okr_query_lst)} okr cnv record!')
    return okr_query_lst


def annotate_sv(vcf, out='fusion.final.txt', other_target_genes=None, okr_fusion_list=None, okr_fusion_pairs=None):
    if other_target_genes:
        targets = set(x.strip() for x in open(other_target_genes))
    else:
        targets = set()

    okr_f_s_d = dict()
    if okr_fusion_list:
        with open(okr_fusion_list) as f:
            for line in f:
                gene, fname = line.strip().split()
                okr_f_s_d[gene] = line.strip()

    okr_f_p_d = dict()
    if okr_fusion_pairs:
        with open(okr_fusion_pairs) as f:
            for line in f:
                gene_pair, fname = line.strip().split()
                okr_f_p_d[gene_pair] = line.strip()

    cmd = f"""bcftools filter -i 'FILTER="PASS" & INFO/PRECISE=1' {vcf} -o {vcf}.filtered """
    cmd2 = f"/nfs2/software/SVAFotate/svafotate.py -i {vcf}.filtered -o {vcf}.filtered.svafotate -ccdg /nfs2/software/SVAFotate/ccdg_sv_afs.bed.gz --gnomad /nfs2/software/SVAFotate/gnomad_sv_afs.bed.gz"
    cmd3 = f"java -Xmx4g -jar /nfs2/software/snpEff/snpEff.jar -v hg19 {vcf}.filtered.svafotate > {vcf}.filtered.svafotate.snpeff"
    cmd4 = f"/nfs2/software/simple_sv_annotation/simple_sv_annotation.py {vcf}.filtered.svafotate.snpeff -o {vcf}.filtered.svafotate.snpeff.simple --known_fusion_pairs /nfs2/software/simple_sv_annotation/fusion_pairs.txt"
    cmd5 = f"bcftools annotate -x 'INFO/ANN' -o {vcf}.filtered.svafotate.snpeff.simplex {vcf}.filtered.svafotate.snpeff.simple"
    cmd6 = f"/nfs2/software/svtools/vcfToBedpe -i {vcf}.filtered.svafotate.snpeff.simplex -o {vcf}.filtered.svafotate.snpeff.simplex.bedpe"
    descs = ['filtering', 'svafotate pop freq', 'snpeff annotate', 'simplify annotation', 'remove ANN', 'vcf2Bed']
    #
    if os.path.exists(f'{vcf}.filtered.svafotate.snpeff.simple'):
        os.remove(f'{vcf}.filtered.svafotate.snpeff.simple')
    if os.path.exists(f'{vcf}.filtered.svafotate.snpeff.simplex'):
        os.remove(f'{vcf}.filtered.svafotate.snpeff.simplex')
    print(f'{vcf}.filtered.svafotate.snpeff.simplex.bedpe')
    if os.path.exists(f'{vcf}.filtered.svafotate.snpeff.simplex.bedpe'):
        # os.remove(f'{vcf}.filtered.svafotate.snpeff.simplex.bedpe')
        print('sv注释结果已经存在，跳过注释')
    else:
        # calculate
        for each, desc in zip([cmd, cmd2, cmd3, cmd4, cmd5, cmd6], descs):
            print('running:', desc)
            check_call(each, shell=True)
    # remove unwanted files
    if os.path.exists(f'{vcf}.filtered'):
        os.remove(f'{vcf}.filtered')
    if os.path.exists(f'{vcf}.filtered.svafotate.snpeff.simple'):
        os.remove(f'{vcf}.filtered.svafotate.snpeff.simple')
    if os.path.exists(f'{vcf}.filtered.svafotate.snpeff'):
        os.remove(f'{vcf}.filtered.svafotate.snpeff')
    if os.path.exists(f'{vcf}.filtered.svafotate'):
        os.remove(f'{vcf}.filtered.svafotate')

    okr_query = []
    with open(f"{vcf}.filtered.svafotate.snpeff.simplex.bedpe") as f, open(out, 'w') as fw:
        #0:CHROM_A 1:START_A 2:END_A 3:CHROM_B 4:START_B 5:END_B 6:ID 7:QUAL 8:STRAND_A 9:STRAND_B 10:TYPE 11:FILTER 12:INFO 13:FORMAT 14:T190079D1L18 15:B180518G1L2
        header = f.readline()
        fw.write('fusion_genes\tbreak_positions\ttype\tOKR_Name\treport\n')
        for line in f:
            lst = line.strip().split('\t')
            if lst[11].lower() not in ['pass']:
                continue
            # SIMPLE_ANN=DEL|EXON_DEL|GNA14|NM_004297.3|Exon1-2del|3,DEL|EXON_DEL|GNA14-AS1|NR_121184.1|Exon2-5del|3;SV_HIGHEST_TIER=3
            # SIMPLE_ANN=INV|GENE_FUSION|ALK&EML4|ALK|KNOWN_FUSION|1
            # SIMPLE_ANN=BND|GENE_FUSION&FRAMESHIFT_VARIANT|ENOX1&TYRO3|NM_006293.3|KNOWN_FUSION|1
            # print(lst[12])
            if 'SIMPLE_ANN=' in lst[12]:
                annot = lst[12].split('SIMPLE_ANN=')[1].split(',')[0].split('|')
                # print(annot)
                if 'GENE_FUSION' not in annot[1] or '&' not in annot[2]:
                    continue
                # print(annot)
                g1, g2 = annot[2].split('&')
                print(g1, g2)
                report = False
                okr_name = ''
                if g1+'-'+g2 in okr_f_p_d:
                    report = True
                    okr_query.append(okr_f_p_d[g1 + '-' + g2])
                    okr_name = okr_query[-1]
                elif g2+'-'+g1 in okr_f_p_d:
                    report = True
                    okr_query.append(okr_f_p_d[g2 + '-' + g1])
                    okr_name = okr_query[-1]
                elif g1 in okr_f_s_d or g2 in okr_f_s_d:
                    report = True
                    if g1 in okr_f_s_d:
                        okr_query.append(okr_f_s_d[g1])
                        okr_name = okr_query[-1]
                    if g2 in okr_f_s_d:
                        okr_query.append(okr_f_s_d[g2])
                        okr_name = okr_query[-1]

                elif g1 in targets or g2 in targets:
                    report = True

                if report:
                    b1 = f'{lst[0]}:{(int(lst[1])+int(lst[2]))//2+1}'
                    b2 = f'{lst[3]}:{(int(lst[4])+int(lst[5]))//2+1}'
                    break_pos = f'{b1}-{b2}'
                    fw.write(f'{g1}--{g2}\t{break_pos}\t{lst[10]}\t{okr_name}\tyes\n')
            else:
                pass
    print(okr_query)
    return okr_query


def filter_germline(vcf, target_genes, comm_trans, genome):
    targets = set(x.strip().split('\t')[0] for x in open(target_genes))
    if not os.path.exists(vcf[:-3]+'txt'):
        print(f'annoated file of {vcf} already exist and skip!')
        annovar_vcf, annovar_txt = annovar_annotation(vcf)
    else:
        annovar_txt = vcf[:-3]+'txt'
    a = pd.read_table(annovar_txt)
    f = ['pathogenic' in x.lower() for x in a['CLNSIG']]
    f0 = [x == 'reviewed_by_expert_panel' for x in a['CLNREVSTAT']]
    if 'Func_refGeneWithVer' in a:
        f1 = [x != 'intronic' for x in a['Func_refGeneWithVer']]
    else:
        f1 = [x != 'intronic' for x in a['Func_refGene']]
    if 'Gene_refGeneWithVer' in a:
        f2 = [x in targets for x in a['Gene_refGeneWithVer']]
    else:
        f2 = [x in targets for x in a['Gene_refGene']]
    # AF > 0.2
    last_col = a.columns[-1]
    # 下面一行试图寻找AF，有可能vcf中没有直接AF信息
    f3 = [any(float(y) >= 0.2 for y in x.split(':')[2].split(',')) >= 0.2 for x in a[last_col]]
    f4 = [x == '.' or float(x) <= 0.01 for x in a['esp6500siv2_all']]
    f5 = [x == '.' or float(x) <= 0.01 for x in a['1000g2015aug_all']]
    f6 = [x == '.' or float(x) <= 0.01 for x in a['ExAC_ALL']]
    f7 = [x == '.' or float(x) <= 0.01 for x in a['gnomAD_exome_ALL']]
    f8 = [x == '.' or float(x) <= 0.01 for x in a['gnomAD_exome_EAS']]
    f9 = [x == '.' or float(x) <= 0.01 for x in a['ExAC_EAS']]
    m = pd.DataFrame([f1, f2, f3, f4, f5, f6, f7, f8, f9]).T
    all_match = m.apply(all, axis=1)
    p = ( pd.Series(f) & pd.Series(f0) & pd.Series(f2) & pd.Series(f3) ) | all_match
    r = a.loc[p]
    if r.shape[0] == 0:
        print(f'{annovar_txt} is empty after filtering!')
    else:
        # 加工一下过滤结果，主要是为了提取common转录本的hgvs信息
        old = annovar_txt[:-3]+'filtered.txt'
        r.to_csv(old)
        out = process_annovar_txt(old, comm_trans, genome=genome, hots=None, af=0.0, bp_num=25)
        os.remove(old)
        return out


def pipeline(input_dir, af=0.02, not_hot_af=0.05, msi_cutoff=10, tmb_cutoff=10,
             hots='/nfs2/database/Panel800/hotspot.xlsx',
             genome='/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta',
             comm_trans='/nfs2/database/Panel800/common_transcripts.txt',
             target_fusion_genes='/nfs2/database/Panel800/other.fusion.gene.list',
             okr_fusion_list='/nfs2/database/Panel800/okr.fusion.gene.list',
             okr_fusion_pairs='/nfs2/database/Panel800/okr.fusion.pair',
             okr_cnv_genes='/nfs2/database/Panel800/okr.cnv.gene.list',
             other_cnv_genes='/nfs2/database/Panel800/other.cnv.gene.list',
             target_germline_genes='/nfs2/database/Panel800/germline.target.txt'
             ):
    """
    通过爬虫的方式注释进行okr注释和爬取报告
    :return:
    """
    # indel processing
    print('Step1: annovar注释'.format(af))
    indel_vcfs = glob(os.path.join(input_dir, 'T*VariantCaller.mutations.*.vcf'))
    if not indel_vcfs:
        raise Exception('Cannot find somatic indel vcf!')
    indel_vcf = sorted(indel_vcfs, key=lambda x:len(x))[0]
    annovar_vcf, annovar_txt = annovar_annotation(indel_vcf)
    tumour_ind = 1 if indel_vcf.startswith('TO') else 0
    _ = filter_vcf_by_af(annovar_vcf, af=af, tumour=tumour_ind)

    print('Step2: 依据{}提取hotspot, 包括af过滤:hotspot af>={} 或 非hotspot af>={}'.format(hots, af, not_hot_af))
    detected = process_annovar_txt(annovar_txt, comm_trans, genome=genome, hots=hots, af=af, not_hot_af=not_hot_af)

    print('Step3: 提取tmb和msi信息并判定是否tmb >= {}, msi >= {}'.format(msi_cutoff, tmb_cutoff))
    tmb_file = glob(os.path.join(input_dir, 'TMB.T*VariantCaller.mutations.*.vcf.txt'))[0]
    with open(tmb_file) as f:
        tmb = float(f.readline().split('TMB:')[1].strip())
    msi_file = glob(os.path.join(input_dir, '*.T*MSI.output'))[0]
    with open(msi_file) as f:
        _ = f.readline()
        msi = float(f.readline().split()[2])

    print('Step4: 提取cnv信息，并判定OKR是否可以注释提取的CNV')
    # cnv processing
    cnv_file = glob(os.path.join(input_dir, '*.HQ20.cnr'))[0]
    out = os.path.join(input_dir, '03.cnv.final.txt')
    cnv_query = parse_cnr(cnv_file, out, okr_targets=okr_cnv_genes, other_target_genes=other_cnv_genes)

    print('Step5: 提取融合信息，并判定OKR是否可以查询该融合信息')
    # fusion processing
    sv_file = glob(os.path.join(input_dir, '*SV.*.vcf'))[0]
    out = os.path.join(input_dir, '04.fusion.final.txt')
    fusion_okr_query = annotate_sv(sv_file, out, target_fusion_genes, okr_fusion_list, okr_fusion_pairs)

    print("Step6: 输出可以用于OKR 查询的列表")
    okr_query_lst = []
    with open('05.TMB_MSI.xls', 'w') as f:
        f.write('metric\tvalue\tOKR_Name\treport\n')
        if msi >= msi_cutoff:
            okr_query_lst.append('Microsatellite instability-High')
            f.write(f'MSI\t{msi}\tMicrosatellite instability-High\tyes\n')
        else:
            f.write(f'MSI\t{msi}\t\tyes\n')
        if tmb >= tmb_cutoff:
            okr_query_lst.append('Tumor Mutational Burden')
            f.write(f'TMB\t{tmb}\tTumor Mutational Burden\tyes\n')
        else:
            f.write(f'TMB\t{tmb}\t\tyes\n')

    if os.path.exists(detected):
        detected = pd.read_excel(detected)
        okr_names = detected['OKR_Name']
        okr_query_lst += [x for x in okr_names if type(x) == str]
    okr_query_lst += cnv_query
    okr_query_lst += fusion_okr_query
    print('okr query:', okr_query_lst)
    with open('06.okr_query.list', 'w') as f:
        for each in okr_query_lst:
            f.write(each + '\n')

    print('Step7: 过滤germline的SNP/indel')
    # germline processing
    germline_vcfs = glob(os.path.join(input_dir, 'germline.*.vcf'))
    if not germline_vcfs:
        print('Cannot find germline indel vcf!')
    else:
        germline_vcf = sorted(germline_vcfs, key=lambda x: len(x))[0]
        filter_germline(germline_vcf, target_germline_genes, comm_trans, genome)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=[
        'pipeline', 'parse_cnr', 'annovar_annotation', 'annotate_sv',
        'filter_germline', 'process_annovar_txt', 'extract_hots',
        'annot_vcf_by_txt', 'filter_vcf_by_af'])

