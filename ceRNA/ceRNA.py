"""
现在形成的circRNA的分析思路:
1.使用STAR比对
2.使用CIRCexplorer2基于比对结果鉴定circRNA
3.根据坐标到数据库CSCD找到其靶向的miRNA，同时也可以找到其可能结合的蛋白。这一步也可以用数据库https.//circinteractome.nia.nih.gov/达到
4.circRNA结合蛋白的资料：https://www.bilibili.com/video/BV1mt4y1a7Gx
5.根据miRNA id 到数据库miRTar获得miRNA 靶向的基因, 并且和差异基因取交集
6.根据上述结果构建互作网络:protein <- circRNA -> miRNA -> gene

思考：
一个circRNA可能结合多个rna结合蛋白
一个circRNA也可以充当海绵靶向多个miRNA
一个miRNA可以靶向上百成千的基因
一个基因可能被多个miRNA靶向
构建完互作网络后，我们的网络可能非常大，如何对可能对circRNA（间接）靶向的基因进行排序？
假设我们发现一个显著下调的circRNA：
1. circRNA靶向的结合蛋白，一般不多
2. miRNA如果仅仅被我们的circRNA靶向，那这个miRNA当然很重要，但通常，miRNA可能和多个circRNA存在互作关系。
所以可以根据miRNA能够结合的circRNA的数目进行排序
3. 同上，一个基因可以被多个miRNA靶向，这时也可以根据一个基因能够被miRNA靶向的基因的数目进行排序
4. 可能还要做一下蛋白互作的网络。


"""


def circRNA2miRNA(circRNAs, cir2mir_db, mir2gene_db=None, out='network.txt',
                  target_gene_candidates=None, geneid2symbol=None):
    """

    :param circRNAs: 文件为4列，chr,start,end,name
    :param cir2mir_db: /nfs2/database/CircRNA/CSCD/hg38_common_circrna_mre.txt
    :param mir2gene_db: /nfs2/database/CircRNA/miRTar/interactions.all
    :param geneid2symbol:
    :param target_gene_candidates: 如果提供，取出mirna靶向的基因与该候选的交集
    :param out:
    :return:
    """
    # 提取待查询circRNA信息
    circ_rnas = dict()
    with open(circRNAs, 'r') as f:
        _header = f.readline()
        for line in f:
            lst = line.strip().split()
            circ_id = lst[0]+':'+lst[1]+'|'+lst[2]
            circ_rnas[circ_id] = lst[3]

    # 提取circRNA能够靶向的miRNA信息
    cir2mir = dict()
    mirna_set = set()
    raw_mir_ids = set()  # 用于后续查询能够靶向某个miRNA的circRNA数目
    with open(cir2mir_db, 'r') as f:
        for line in f:
            lst = line.split()
            circ = lst[0]
            mir = lst[1]
            if circ in circ_rnas:
                raw_mir_ids.add(mir)
                circ = circ_rnas[circ]
                if '/' in mir:
                    prefix, ops = mir.split('-', 1)
                    mir_ids = {'hsa-'+prefix+'-'+x for x in ops.split('/')}
                else:
                    mir_ids = {'hsa-'+lst[1]}
                cir2mir.setdefault(circ, set()).update(mir_ids)
                mirna_set.update(mir_ids)

    # 提取miRNA能够靶向的基因信息
    mir2target = dict()
    target_genes = set()
    if mir2gene_db:
        geneid2symbol = dict(x.strip().split()[:2] for x in open(geneid2symbol))
        if target_gene_candidates:
            degs = set(x.strip().split('\t')[0] for x in open(target_gene_candidates))
        else:
            degs = set()
        with open(mir2gene_db, 'r') as f:
            for line in f:
                lst = line.split()
                if lst[0] in mirna_set:
                    genes = lst[1].split(";")
                    if degs:
                        genes = set(genes) & degs
                        target_genes.update(genes)
                    if genes:
                        gene_symbols = [geneid2symbol[x] if x in geneid2symbol else x for x in sorted(genes)]
                        mir2target[lst[0]] = gene_symbols

    # # 从数据库中提取mirna能够靶向的circRNA数目信息
    # mir2cir_num = dict()
    # with open(cir2mir_db, 'r') as f:
    #     for line in f:
    #         lst = line.split()
    #         circ = lst[0]
    #         mir = lst[1]
    #         if mir in raw_mir_ids:
    #             if '/' in mir:
    #                 prefix, ops = mir.split('-', 1)
    #                 mir_ids = {'hsa-'+prefix+'-'+x for x in ops.split('/')}
    #             else:
    #                 mir_ids = {'hsa-'+lst[1]}
    #             for each in mir_ids:
    #                 mir2cir_num.setdefault(each, 0)
    #                 mir2cir_num[each] += 1
    # mir2cir_num = dict(sorted(zip(mir2cir_num.keys(), mir2cir_num.values()), key=lambda x:x[1]))
    # with open('mir2circRNA_num.txt', 'w') as f:
    #     for k, v in mir2cir_num.items():
    #         f.write(f'{k}\t{v}\n')

    # # 从数据库中提取基因被靶向的miRNA数量
    # gene2mir_num = dict()
    # with open(mir2gene_db, 'r') as f:
    #     for line in f:
    #         lst = line.split()
    #         genes = lst[1].split(";")
    #         genes = set(genes) & target_genes
    #         for gene in genes:
    #             gene = geneid2symbol[gene]
    #             gene2mir_num.setdefault(gene, 0)
    #             gene2mir_num[gene] += 1
    # gene2mir_num = dict(sorted(zip(gene2mir_num.keys(), gene2mir_num.values()), key=lambda x:x[1]))
    # with open('gene2miRNA_num.txt', 'w') as f:
    #     for k, v in gene2mir_num.items():
    #         f.write(f'{k}\t{v}\n')

    # 输出互作网络，网络的节点属性需另外准备
    with open(out, 'w') as f:
        for circ, mirnas in cir2mir.items():
            for mir in mirnas:
                f.write(f'{circ}\t{mir}\n')

        for mir, genes in mir2target.items():
            for gene in genes:
                f.write(f'{mir}\t{gene}\n')

    if target_genes:
        with open('target_gene.list', 'w') as f:
            _ = [f.write(x+'\n') for x in target_genes]


def mir2gene_TargetScanHuman(cicrRNAs, mir2gene_db, out='mir2gene.network.txt'):
    mirna_set = set()
    with open(cicrRNAs) as f:
        for line in f:
            mir = line.strip().split()[0]
            mirna_set.add(mir)
    print(mirna_set)

    with open(mir2gene_db, 'r') as f, open(out, 'w') as fw:
        for line in f:
            lst = line.split()
            mir = lst[0]
            if '/' in mir:
                prefix, ops = mir.split('-', 1)
                mir_ids = {'hsa-' + prefix + '-' + x for x in ops.split('/')}
            else:
                mir_ids = {'hsa-' + mir}
            # print(mir_ids)
            for mir in mir_ids:
                if mir in mirna_set:
                    gene_id = lst[1].rsplit('.', 1)[0]
                    gene_symbol = lst[2]
                    fw.write(f'{mir}\t{gene_id}\t{gene_symbol}\n')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
