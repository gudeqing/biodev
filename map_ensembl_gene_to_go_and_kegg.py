# coding=utf-8
import os
import wget
import gzip
import logging

"""
目的：获取某物种的ensembl_id对应的go注释和kegg注释信息，需要联网获取
根据ftp://ftp.ncbi.nih.gov/gene/DATA/，该网址更新信息频繁，最好每次都执行该脚本。
 ensembl id --> ncbi_id
 ncbi_id --> go_id
 ncbi_id --> kegg k_id
 ncbi_id --> kegg path
 ncbi_id --> kegg enzyme
 Note that KEGG IDs are the same as Entrez Gene IDs for most species anyway
 KEGG uses Entrez Gene ID as its standard gene ID
注释：经检查发现，至少在human中，kegg中特定物种的基因id就是ncbi_id前加上物种缩写名，如'1'对应的kegg中的"hsa:1"
"""


def set_logger(name='log.info', logger_id='x'):
    logger = logging.getLogger(logger_id)
    fh = logging.FileHandler(name, mode='w+')
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    return logger


def ensembl2ncbi(gene2ensembl_file="gene2ensembl.gz", tax_id='9606'):
    """
    #tax_id	GeneID	Ensembl_gene_identifier	RNA_nucleotide_accession.version	Ensembl_rna_identifier	protein_accession.version	Ensembl_protein_identifier
    7227	30970	FBgn0040373	NM_001297803.1	FBtr0340207	NP_001284732.1	FBpp0309182
    7227	30970	FBgn0040373	NM_130477.4	FBtr0070108	NP_569833.1	FBpp0070103
    7227	30970	FBgn0040373	NM_166834.2	FBtr0070107	NP_726658.1	FBpp0070102
    7227	30971	FBgn0040372	NM_001272159.1	FBtr0332992	NP_001259088.1	FBpp0305207
    """
    # human taxonomy id is 9606
    logger = set_logger()
    if not os.path.exists(gene2ensembl_file):
        print('downloading gene2ensembl.gz from ncbi')
        wget.download('ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz', gene2ensembl_file, bar=None)
    gene2ncbi = dict()
    transcript2ncbi = dict()
    out_name = '{}.ensembl.gene2ncbi.txt'.format(tax_id if tax_id != '9606' else 'hsa')
    out_name2 = '{}.ensembl.trans2ncbi.txt'.format(tax_id if tax_id != '9606' else 'hsa')
    f = open(out_name, 'w')
    f2 = open(out_name2, 'w')
    for line in gzip.open(gene2ensembl_file, mode='rt', encoding='utf-8'):
        if line.startswith(tax_id):
            line_list = line.strip().split('\t')
            if line_list[2] in gene2ncbi:
                if line_list[1] not in gene2ncbi[line_list[2]]:
                    if line_list[2].strip() == '-':
                        continue
                    logger.info(line_list[2]+' has multiple ncbi gene id')
                    gene2ncbi[line_list[2]].append(line_list[1])
                    f.write('{}\t{}\n'.format(line_list[2], line_list[1]))
            else:
                gene2ncbi[line_list[2]] = [line_list[1]]
                f.write('{}\t{}\n'.format(line_list[2], line_list[1]))

            if line_list[4] in transcript2ncbi:
                if line_list[1] not in transcript2ncbi[line_list[4]]:
                    if line_list[4].strip() == '-':
                        continue
                    logger.info(line_list[4]+' has multiple ncbi gene id')
                    transcript2ncbi[line_list[4]].append(line_list[1])
                    f2.write('{}\t{}\n'.format(line_list[4], line_list[1]))
            else:
                transcript2ncbi[line_list[4]] = [line_list[1]]
                f2.write('{}\t{}\n'.format(line_list[4], line_list[1]))
    f.close()
    f2.close()
    return gene2ncbi, transcript2ncbi


def ncbi2go(gene2go_file='gene2go.gz', tax_id='9606'):
    if not os.path.exists(gene2go_file):
        print('downloading gene2go.gz from ncbi')
        wget.download('ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz', gene2go_file, bar=None)
    gene2go = dict()
    out_name = '{}.gene2go.txt'.format(tax_id if tax_id != '9606' else 'hsa')
    f = open(out_name, 'w')
    for line in gzip.open(gene2go_file, mode='rt', encoding='utf-8'):
        if line.startswith(tax_id):
            line_list = line.split()
            gene2go.setdefault(line_list[1], list())
            gene2go[line_list[1]].append(line_list[2])
            f.write('{}\t{}\n'.format(line_list[1], line_list[2]))
    f.close()
    return gene2go


def ncbi2ko(species='hsa'):
    ncbi2kegg_url = 'http://rest.kegg.jp/link/ko/{}'.format(species)
    ncbi2kegg_file = '{}.gene2kegg.txt'.format(species)
    print('downloading {} ko information from kegg'.format(species))
    wget.download(ncbi2kegg_url, ncbi2kegg_file, bar=None)
    n2k = dict()
    with open(ncbi2kegg_file) as f:
        for line in f:
            ncbi_id, ko_id = line.strip().split()
            ncbi_id = ncbi_id.split(":")[1]
            ko_id = ko_id.split(":")[1]
            n2k.setdefault(ncbi_id, list())
            n2k[ncbi_id].append(ko_id)
    return n2k


def ncbi2path(species='hsa'):
    url = 'http://rest.kegg.jp/link/pathway/{}'.format(species)
    print('downloading {} pathway information from kegg'.format(species))
    g2p_file = '{}.gene2path.txt'.format(species)
    wget.download(url, g2p_file, bar=None)
    g2p = dict()
    with open(g2p_file) as f:
        for line in f:
            gene, path = line.strip().split()
            gene = gene.strip().split(":", 1)[1]  # 把gene前的物种缩写去掉，得到的就是ncbi，至少hsa是这样的
            path = path.strip()
            g2p.setdefault(gene, list())
            g2p[gene].append(path)
    return g2p


def ncbi2enzyme(species='hsa'):
    url = 'http://rest.kegg.jp/link/enzyme/{}'.format(species)
    g2e_file = '{}.gene2enzyme.txt'.format(species)
    print('downloading {} enzyme information from kegg'.format(species))
    wget.download(url, g2e_file, bar=None)
    g2e = dict()
    with open(g2e_file) as f:
        for line in f:
            gene, enzyme = line.strip().split()
            gene = gene.strip().split(":", 1)[1]  # 把gene前的物种缩写去掉，得到的就是ncbi，至少hsa是这样的
            enzyme = enzyme.strip()
            g2e.setdefault(gene, list())
            g2e[gene].append(enzyme)
    return g2e


def ko2path():
    _ = 'not used now'
    k2p_url = 'http://rest.kegg.jp/link/pathway/ko'
    k2p_file = 'ko2pathway.txt'
    wget.download(k2p_url, k2p_file, bar=None)
    k2p_dict = dict()
    with open(k2p_file) as f:
        for line in f:
            ko, path = line.strip().split()
            if path.startswith('path:map'):
                continue
            ko = ko.split(':')[1]
            path = path.split(':')[1]
            k2p_dict.setdefault(ko, list())
            k2p_dict[ko].append(path)
    return k2p_dict


def map_ensembl_to_go_kegg(species='hsa', tax_id='9606'):
    g2n, t2n = ensembl2ncbi(tax_id=tax_id)
    n2go = ncbi2go(tax_id=tax_id)
    n2k = ncbi2ko(species=species)
    n2p = ncbi2path(species=species)
    n2e = ncbi2enzyme(species=species)
    # map gene info
    gene2go = dict()
    gene2path = dict()
    gene2ko = dict()
    gene2enzyme = dict()
    for g, n_lst in g2n.items():
        for n in n_lst:
            if n in n2go:
                gene2go.setdefault(g, list())
                gene2go[g] += n2go[n]
            if n in n2p:
                gene2path.setdefault(g, list())
                gene2path[g] += n2p[n]
            if n in n2k:
                gene2ko.setdefault(g, list())
                gene2ko[g] += n2k[n]
            if n in n2e:
                gene2enzyme.setdefault(g, list())
                gene2enzyme[g] += n2e[n]
    # map transcript info
    trans2go = dict()
    trans2path = dict()
    trans2ko = dict()
    trans2enzyme = dict()
    for g, n_lst in t2n.items():
        for n in n_lst:
            if n in n2go:
                trans2go.setdefault(g, list())
                trans2go[g] += n2go[n]
            if n in n2p:
                trans2path.setdefault(g, list())
                trans2path[g] += n2p[n]
            if n in n2k:
                trans2ko.setdefault(g, list())
                trans2ko[g] += n2k[n]
            if n in n2e:
                trans2enzyme.setdefault(g, list())
                trans2enzyme[g] += n2e[n]

    result = zip(
        [species+'.ensembl.'+x for x in [
            'gene2go.txt', 'gene2ko.txt', 'gene2path.txt', 'gene2enzyme.txt',
            'trans2go.txt', 'trans2ko.txt', 'trans2path.txt', 'trans2enzyme.txt']],
        [gene2go, gene2ko, gene2path, gene2enzyme,
         trans2go, trans2ko, trans2path, trans2enzyme]
    )
    for name, data in result:
        with open(name, 'w') as f:
            for k, v in data.items():
                for each_v in v:
                    f.write('{}\t{}\n'.format(k, each_v))
    return


if __name__ == '__main__':
    map_ensembl_to_go_kegg(species='hsa', tax_id='9606')
    # wget.download('http://rest.kegg.jp/link/enzyme/ko', 'ko2enzyme.txt')
    from urllib.request import urlretrieve
    url = 'https://www.kegg.jp/kegg-bin/download_htext?htext=br08901&format=htext&filedir='
    print('downloading br08901.keg')
    urlretrieve(url, 'br08901.keg')
    if not os.path.exists('go-basic.obo'):
        print('downloading go-basic.obo')
        wget.download('http://purl.obolibrary.org/obo/go/go-basic.obo', 'go-basic.obo', bar=None)

