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


def parse_class_mutation(infile):
    table = pd.read_csv(infile, header=0, index_col=0, sep=None, engine='python')
    pattern = re.compile(r'([A-Z0-9]+) mutation$')
    genes = set()
    for each in table['variant_class']:
        match = pattern.match(each)
        if match:
            genes.add(match.groups()[0])
    return genes


def parse_class_aberration(infile):
    table = pd.read_csv(infile, header=0, index_col=0, sep=None, engine='python')
    pattern = re.compile(r'([A-Z0-9]+) aberration$')
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


def max_distance_to_root(tree, leaf):
    path = []

    def find_root(tree, leaf):
        """"
        tree = dict(a='b', b=['c', 'f'], c=['d', 'g'], d=['e', 'r'])
                  f   e
                 /   /
        a - b - c - d
                 \   \
                  g   r
        """
        nonlocal path
        root = leaf
        if root in tree:
            for root in tree[root]:
                path.append(root)
                yield from find_root(tree, root)
        else:
            # 除了已经返回过的节点，每次返回是最深的节点，即距离当前leaf最远的节点
            yield root, len(path), path
    return next(find_root(tree, leaf))[1]


def mimic_okr(hot_table, class_condition, parent_info, out='okr_mutation.xlsx'):
    class2mut, mut2class = parse_class_conidition(class_condition)
    print(f'有{len(mut2class)}具体突变被分为{len(class2mut)}类')
    exon2class = parse_class_exon(class_condition)
    protein2class = parse_class_half_known(class_condition)
    deleterious_genes = parse_class_deleterious(class_condition)
    activating_genes = parse_class_activating(class_condition)
    oncokb_genes = oncokb_truncating_genes()
    splice_genes = parse_class_splice(class_condition)
    mutation_genes = parse_class_mutation(class_condition)
    aberration_genes = parse_class_aberration(class_condition)
    print('deleterious_genes:{}, activating_genes:{}, splice_genes:{}, mutation_genes:{}, aberration_genes:{}'.format(
        len(deleterious_genes), len(activating_genes), len(splice_genes), len(mutation_genes), len(aberration_genes)))

    table = pd.read_csv(hot_table, header=0, sep=None, engine='python').fillna('.')
    table = table.set_index(['#CHROM', 'POS', 'REF', 'ALT'])
    targets = ['Gene_refGene', 'ExonNumber', 'Final_pHGVS', 'Variation Type',
               'Oncogenicity', 'Func_refGene', 'ExonicFunc_refGene', 'CLNREVSTAT', 'CLNSIG',]
    target_df = table.loc[:, targets]
    half = re.compile(r'p.([A-Z][0-9]+[*]?).*')
    result = list()
    for ind, (gene, exon, phgvs, tp, onco, func, exonic_func, clnrevstat, clnsig) in target_df.iterrows():
        result.append([set(), ''])
        go_on = True

        if go_on:
            # 搜索是否在具体的分类里
            key = gene + ':' + phgvs
            if key in mut2class:
                result[-1][0].update(mut2class[key])
                result[-1][1] = '来自OKR具体分类'
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
                                result[-1][0].add(each)
                                result[-1][1] = '根据半已知突变和数据库注释共同推断分类'
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

                match = onco.lower() in ['Likely_Oncogenic', 'Oncogenic']
                match2 = clnrevstat.lower() in ['reviewed_by_expert_panel'] and ('pathogenic' in clnsig.lower())
                match3 = gene in oncokb_genes
                match4 = 'splicing' in func.lower()
                consider = ['frameshift deletion', 'frameshift insertion',
                            'frameshift substitution', 'stopgain', 'stoploss']
                match5 = exonic_func.strip().lower() in consider
                if match or match2 or (match3 and (match4 or match5)):
                    candidates.append(f'{gene} exon {exon} activating mutation')

                candidates.append(f'{gene} exon {exon} mutation')
                candidates.append(f'{gene} {half.groups()[0]} mutation')
                if set(candidates) & cls_set:
                    for each in candidates:
                        if each in cls_set:
                            result[-1][0].add(each)
                            result[-1][1] = '根据外显子号和突变类型和数据库注释共同推断分类'
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
                    result[-1][0].add(f'{gene} deleterious mutation')
                    result[-1][1] = '根据oncogenic和clinvar判断是否致病'
                elif gene in activating_genes:
                    result[-1][0].add(f'{gene} activating mutation')
                    result[-1][1] = '根据oncogenic和clinvar判断是否致病'
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
                    result[-1][0].add(f'{gene} deleterious mutation')
                    result[-1][1] = '根据oncokb判断是否致病'
                elif gene in activating_genes:
                    result[-1][0].add(f'{gene} activating mutation')
                    result[-1][1] = '根据oncokb判断是否致病'
                else:
                    go_on = True
            else:
                go_on = True

        if go_on:
            # splicing 判断
            go_on = False
            if 'splicing' in func.lower() and (gene in splice_genes):
                result[-1][0].add(f'{gene} splice site mutation')
                result[-1][1] = '根据func注释是否包含splice推断分类'
            else:
                go_on = True

        if go_on:
            # 归为mutation
            go_on = False
            if gene in mutation_genes:
                result[-1][0].add(f'{gene} mutation')
                result[-1][1] = '根据OKR mutation'
            else:
                go_on = True

        if go_on:
            # 归为mutation
            go_on = False
            if gene in aberration_genes:
                result[-1][0].add(f'{gene} aberration')
                result[-1][1] = '根据OKR aberration'
            else:
                go_on = True

    # 对于多个分类的情况，根据parent信息筛选出层级最低的那个
    tree = parse_class_parent(parent_info)
    all_cls = {x.strip().split('\t')[1] for x in open(parent_info)}
    cls_lst = []
    reason_lst = []
    for cls_set, reason in result:
        if len(cls_set) > 1:
            cls_depth = map(max_distance_to_root, [tree for x in range(len(cls_set))], list(cls_set))
            selected = sorted(zip(list(cls_set), cls_depth), key=lambda x:x[1])[-1][0]

            print(selected, 'vs ', cls_set)
            cls_lst.append(selected)
        else:
            if cls_set:
                for each in cls_set:
                    if each not in all_cls:
                        raise Exception(f'分类名{each}不在已知分类里!')
                    cls_lst.append(each)
            else:
                cls_lst.append('')
        reason_lst.append(reason)

    target_df['OKR_Name'] = cls_lst
    target_df['myReason'] = reason_lst
    target_df.to_excel(out, merge_cells=False)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
