import os
import re
import pandas as pd


def parse_class_conidition(infile):
    table = pd.read_csv(infile, header=0, index_col=0, sep=None, engine='python')
    class2mutation = dict()
    mutation2class = dict()
    for ind, row in table.iterrows():
        variant_class = row['variant_class']
        c_lst = eval(row['condition'])
        genes = []
        proteins = []
        for c in c_lst:
            if 'field' in c and c['field'] == 'gene':
                gene = c['value']
                if ',' in gene:
                    genes += gene.split(',')
                else:
                    genes.append(gene)
            elif 'field' in c and c['field'] == 'protein':
                if c['operator'] == 'in':
                    proteins += c['value'].split(',')
                else:
                    proteins.append(c['value'])
        if genes and proteins:
            class2mutation.setdefault(variant_class, set())
            for gene in genes:
                for protein in proteins:
                    key = gene+':'+protein
                    class2mutation[variant_class].add(key)
                    mutation2class.setdefault(key, set())
                    mutation2class[key].add(variant_class)
    return class2mutation, mutation2class


def parse_class_exon(infile):
    table = pd.read_csv(infile, header=0, index_col=0, sep=None, engine='python')
    # EGFR exon 19 activating mutation
    possible_mode = [
        "deletion",
        "insertion",
        "resistance mutation",
        "sensitizing mutation",
        "activating mutation",
        "skipping",
        "mutation",
    ]
    exon_dict = dict()
    pattern = re.compile(r'(.*) exon ([0-9]+) (.*)')
    for each in table['variant_class']:
        match = pattern.match(each)
        if match:
            gene, exon_num, mode = match.groups()
            key = gene + ':' + exon_num
            exon_dict.setdefault(key, set())
            exon_dict[key].add(each)
    return exon_dict


def parse_class_half_known(infile):
    table = pd.read_csv(infile, header=0, index_col=0, sep=None, engine='python')
    possible_mode = [
        "deletion",
        "frame shift",
        "mutation",
        "mutation status", # "oncomine_variant_class" -> "Hotspot,Deleterious"}
    ]
    half_known = dict()
    pattern = re.compile(r'(.*) ([A-Z][0-9]+[*]?) (.*)')
    for each in table['variant_class']:
        match = pattern.match(each)
        if match:
            try:
                gene, half_mutation, mode = match.groups()
            except Exception as e:
                print(each)
                exit(e)
            key = gene + ':' + half_mutation
            half_known.setdefault(key, set())
            half_known[key].add(each)
    return half_known


def parse_class_deleterious(infile):
    table = pd.read_csv(infile, header=0, index_col=0, sep=None, engine='python')
    pattern = re.compile(r'([A-Z0-9]+) deleterious mutation')
    genes = set()
    for each in table['variant_class']:
        match = pattern.match(each)
        if match:
            genes.add(match.groups()[0])
    return genes


def parse_class_activating(infile):
    table = pd.read_csv(infile, header=0, index_col=0, sep=None, engine='python')
    pattern = re.compile(r'([A-Z0-9]+) activating mutation')
    genes = set()
    for each in table['variant_class']:
        match = pattern.match(each)
        if match:
            genes.add(match.groups()[0])
    return genes


def parse_class_splice(infile):
    table = pd.read_csv(infile, header=0, index_col=0, sep=None, engine='python')
    pattern = re.compile(r'([A-Z0-9]+) splice site mutation')
    genes = set()
    for each in table['variant_class']:
        match = pattern.match(each)
        if match:
            genes.add(match.groups()[0])
    return genes


def oncokb_truncating_genes():
    genes = "ACTG1 AMER1 ANKRD11 APC ARID1A ARID1B ARID2 ARID3A ARID4A ARID4B ARID5B ASXL1 ASXL2 ATM ATP6V1B2 ATR ATRX ATXN2 AXIN1 AXIN2 B2M BACH2 BAP1 BARD1 BBC3 BCL10 BCL11B BCL2L11 BCOR BCORL1 BIRC3 BLM BMPR1A BRCA1 BRCA2 BRIP1 BTG1 CASP8 CBFB CBL CD58 CDC73 CDH1 CDK12 CDKN1A CDKN1B CDKN2A CDKN2B CDKN2C CEBPA CHEK1 CHEK2 CIC CIITA CRBN CREBBP CTCF CUX1 CYLD DAXX DDX3X DICER1 DNMT3A DNMT3B DTX1 DUSP22 DUSP4 ECT2L EED EGR1 ELF3 EP300 EP400 EPCAM EPHA3 EPHA7 ERCC2 ERCC3 ERCC4 ERF ERRFI1 ESCO2 EZH2 FAM175A FAM58A FANCA FANCC FANCD2 FAS FAT1 FBXO11 FBXW7 FH FLCN FOXA1 FOXL2 FOXO1 FOXP1 FUBP1 GATA3 GPS2 GRIN2A HIST1H1B HIST1H1D HLA-A HLA-B HLA-C HNF1A HOXB13 ID3 INHA INPP4B INPPL1 IRF1 IRF8 JAK1 JARID2 KDM5C KDM6A KEAP1 KLF2 KLF4 KMT2A KMT2B KMT2C KMT2D LATS1 LATS2 LZTR1 MAP2K4 MAP3K1 MAX MED12 MEN1 MGA MITF MLH1 MRE11A MSH2 MSH3 MSH6 MUTYH NBN NCOR1 NF1 NF2 NFE2 NFKBIA NKX3-1 NOTCH1 NOTCH4 NPM1 NSD1 NTHL1 PALB2 PARK2 PARP1 PAX5 PBRM1 PHF6 PHOX2B PIK3R1 PIK3R2 PIK3R3 PMAIP1 PMS1 PMS2 POT1 PPP2R1A PPP6C PRDM1 PTCH1 PTEN PTPN1 PTPN2 PTPRD PTPRS PTPRT RAD21 RAD50 RAD51 RAD51B RAD51C RAD51D RASA1 RB1 RBM10 RECQL RECQL4 REST RNF43 RTEL1 RUNX1 RYBP SDHA SDHAF2 SDHB SDHC SDHD SESN1 SESN2 SESN3 SETD2 SH2B3 SH2D1A SHQ1 SLX4 SMAD2 SMAD3 SMAD4 SMARCA4 SMARCB1 SOCS1 SOX17 SOX9 SPEN SPOP SPRED1 STAG2 STK11 SUFU SUZ12 TBX3 TCF3 TCF7L2 TET1 TET2 TGFBR1 TGFBR2 TMEM127 TNFAIP3 TNFRSF14 TOP1 TP53 TP53BP1 TP63 TSC1 TSC2 VHL WT1 XRCC2"
    return genes.split()


def parse_class_parent(infile):
    result = dict()
    with open(infile) as f:
        header = f.readline()
        for line in f:
            lst = line.strip().split('\t')
            if len(lst) == 3:
                _, variant_class, parent_variant_class = lst
                result.setdefault(variant_class, set())
                result[variant_class].add(parent_variant_class)
    return result


def relation_graph(infile, outfile='relation.svg', rankdir='LR', class_condition=None):
    import pygraphviz as pgv
    graph = pgv.AGraph(directed=True, rankdir=rankdir)
    class2mut = dict()
    if class_condition:
        class2mut, _ = parse_class_conidition(class_condition)
    for source, targets in parse_class_parent(infile).items():
        for target in targets:
            graph.add_edge(source, target)
            if target in class2mut:
                name = '\n'.join(class2mut[target])
                graph.add_edge(name, target)
                graph.add_node(name, shape='box', color='red')
            if source in class2mut:
                name = '\n'.join(class2mut[source])
                graph.add_edge(name, source)
                graph.add_node(name, shape='box', color='red')
    img_fmt = os.path.splitext(outfile)[1][1:]
    graph.draw(path=outfile, format=img_fmt, prog='dot')


def mimic_okr(hot_table, class_condition, out='okr_mutation.xls'):
    class2mut, mut2class = parse_class_conidition(class_condition)
    exon2class = parse_class_exon(class_condition)
    protein2class = parse_class_half_known(class_condition)
    deleterious_genes = parse_class_deleterious(class_condition)
    activating_genes = parse_class_activating(class_condition)
    oncokb_genes = oncokb_truncating_genes()
    splice_genes = parse_class_splice(class_condition)

    table = pd.read_csv(hot_table, header=0, sep=None, engine='python').fillna('.')
    table = table.set_index(['#CHROM', 'POS', 'REF', 'ALT'])
    targets = ['Gene_refGene', 'ExonNumber', 'Final_pHGVS', 'Variation Type',
               'Oncogenicity', 'Func_refGene', 'ExonicFunc_refGene', 'CLNREVSTAT', 'CLNSIG',]
    target_df = table.loc[:, targets]
    half = re.compile(r'p.([A-Z][0-9]+[*]?).*')
    result = list()
    for ind, (gene, exon, phgvs, tp, onco, func, exonic_func, clnrevstat, clnsig) in target_df.iterrows():
        result.append(set())
        go_on = True

        if go_on:
            # 搜索是否在具体的分类里
            key = gene + ':' + phgvs
            if key in mut2class:
                result[-1].update(mut2class[key])
                go_on = False
            else:
                go_on = True

        if go_on:
            # 搜索是否在半已知突变的分类里
            go_on = False
            match = half.match(phgvs)
            if not match:
                print(f'{ind} {gene} phgvs: {phgvs} cannot be parsed!')
                go_on = True
            else:
                key = gene + ':' + match.groups()[0]
                if key in protein2class:
                    cls_set = protein2class[key]
                    candidates = []
                    if tp.lower() in ['deletion']:
                        candidates.append(f'{gene} {match.groups()[0]} deletion')
                    if exonic_func.startswith('frameshift'):
                        candidates.append(f'{gene} {match.groups()[0]} frame shift')
                    candidates.append(f'{gene} {match.groups()[0]} mutation')
                    if set(candidates) & cls_set:
                        for each in candidates:
                            if each in cls_set:
                                result[-1].add(each)
                                break
                    else:
                        print(f'{ind} {gene} {phgvs} not in range of {cls_set}')
                        go_on = True
                else:
                    go_on = True

        if go_on:
            # 根据外显子号和突变类型分类
            go_on = False
            key = gene + ':' + str(exon)
            if key in exon2class:
                cls_set = exon2class[key]
                candidates = []
                if tp.lower() in ['deletion', 'insertion']:
                    candidates.append(f'{gene} exon {exon} {tp.lower()}')
                if onco.lower() in ['Likely_Oncogenic', 'Oncogenic']:
                    candidates.append(f'{gene} exon {exon} activating mutation')
                candidates.append(f'{gene} exon {exon} mutation')
                candidates.append(f'{gene} {half.groups()[0]} mutation')
                if set(candidates) & cls_set:
                    for each in candidates:
                        if each in cls_set:
                            result[-1].add(each)
                            break
                else:
                    print(f'{ind} {gene} {phgvs} not in range of {cls_set}')
                    go_on = True
            else:
                go_on = True

        if go_on:
            # 根据oncogenic和clinvar判断是否致病
            go_on = False
            match = onco.lower() in ['likely_oncogenic', 'oncogenic']
            match2 = clnrevstat.lower() in ['reviewed_by_expert_panel'] and ('pathogenic' in clnsig.lower())
            if match or match2:
                if gene in deleterious_genes:
                    result[-1].add(f'{gene} deleterious mutation')
                elif gene in activating_genes:
                    result[-1].add(f'{gene} activating mutation')
                else:
                    go_on = True
            else:
                go_on = True

        if go_on:
            # 根据oncokb判断是否致病
            go_on = False
            match = gene in oncokb_genes
            match2 = 'splicing' in func.lower()
            consider = ['frameshift deletion', 'frameshift insertion', 'frameshift substitution', 'stopgain', 'stoploss']
            match3 = exonic_func.strip().lower() in consider
            if match and (match2 or match3):
                if gene in deleterious_genes:
                    result[-1].add(f'{gene} deleterious mutation')
                elif gene in activating_genes:
                    result[-1].add(f'{gene} activating mutation')
                else:
                    go_on = True
            else:
                go_on = True

        if go_on:
            # splicing 判断
            go_on = False
            if 'splicing' in func.lower() and (gene in splice_genes):
                result[-1].add(f'{gene} splice site mutation')
            else:
                go_on = True

    target_df['OKR_Name'] = [";".join(x) for x in result]
    target_df.to_csv(out, sep='\t')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
