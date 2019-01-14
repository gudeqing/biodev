# coding=utf-8
import os
import wget
import gzip

"""
目的：获取某物种的ensembl_id对应的go注释和kegg注释信息，需要联网获取
根据ftp://ftp.ncbi.nih.gov/gene/DATA/，该网址更新信息频繁，最好每次都执行该脚本。
 把ensembl id转换为对应的ncbi_id
 ncbi_id --> go_id
 ncbi_id --> kegg的ko_id
 ko_id --> path_id
注释：经检查发现，至少在human中，kegg中特定物种的基因id就是ncbi_id前加上物种缩写名，如'1'对应的kegg中的"hsa:1"
"""


def ensembl2ncbi(gene2ensembl_file="gene2ensembl.gz", tax_id='9606'):
    """
    #tax_id	GeneID	Ensembl_gene_identifier	RNA_nucleotide_accession.version	Ensembl_rna_identifier	protein_accession.version	Ensembl_protein_identifier
    7227	30970	FBgn0040373	NM_001297803.1	FBtr0340207	NP_001284732.1	FBpp0309182
    7227	30970	FBgn0040373	NM_130477.4	FBtr0070108	NP_569833.1	FBpp0070103
    7227	30970	FBgn0040373	NM_166834.2	FBtr0070107	NP_726658.1	FBpp0070102
    7227	30971	FBgn0040372	NM_001272159.1	FBtr0332992	NP_001259088.1	FBpp0305207
    """
    # human taxonomy id is 9606
    if not os.path.exists(gene2ensembl_file):
        wget.download('ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz', gene2ensembl_file)
    gene2ncbi = dict()
    transcript2ncbi = dict()
    for line in gzip.open(gene2ensembl_file, mode='rt'):
        if line.startswith(tax_id):
            line_list = line.split()
            if line_list[2] in gene2ncbi:
                if line_list[1] not in gene2ncbi[line_list[2]]:
                    print(line_list[2]+' has multiple ncbi gene id')
                    gene2ncbi[line_list[2]].append(line_list[1])
            else:
                gene2ncbi[line_list[2]] = [line_list[1]]
            if line_list[4] in transcript2ncbi:
                if line_list[1] not in transcript2ncbi[line_list[4]]:
                    print(line_list[4]+' has multiple ncbi gene id')
                    transcript2ncbi[line_list[4]].append(line_list[1])
            else:
                transcript2ncbi[line_list[4]] = [line_list[1]]
    return gene2ncbi, transcript2ncbi


def ncbi2go(gene2go_file='gene2go.gz', tax_id='9606'):
    if not os.path.exists(gene2go_file):
        wget.download('ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz', gene2go_file)
    gene2go = dict()
    for line in gzip.open(gene2go_file, mode='rt'):
        if line.startswith(tax_id):
            line_list = line.split()
            gene2go.setdefault(line_list[1], list())
            gene2go[line_list[1]].append(line_list[2])
    return gene2go


def ncbi2ko(species='hsa'):
    ncbi2kegg_url = 'http://rest.kegg.jp/link/ko/{}'.format(species)
    ncbi2kegg_file = 'ncbi2ko.txt'
    wget.download(ncbi2kegg_url, ncbi2kegg_file)
    n2k = dict()
    with open(ncbi2kegg_file) as f:
        for line in f:
            ncbi_id, ko_id = line.strip().split()
            ncbi_id = ncbi_id.split(":")[1]
            ko_id = ko_id.split(":")[1]
            n2k.setdefault(ncbi_id, list())
            n2k[ncbi_id].append(ko_id)
    return n2k


def ko2path():
    k2p_url = 'http://rest.kegg.jp/link/pathway/ko'
    k2p_file = 'ko2pathway.txt'
    wget.download(k2p_url, k2p_file)
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
    k2p = ko2path()
    gene2go = dict()
    gene2path = dict()
    gene2ko = dict()
    for g, n_lst in g2n.items():
        for n in n_lst:
            if n in n2go:
                gene2go.setdefault(g, list())
                gene2go[g] += n2go[n]
            if n in n2k:
                path_lst = list()
                k_lst = n2k[n]
                gene2ko.setdefault(g, list())
                gene2ko[g] += k_lst
                for k in k_lst:
                    if k in k2p:
                        path_lst.extend(k2p[k])
                gene2path.setdefault(g, list())
                gene2path[g] += path_lst

    trans2go = dict()
    trans2path = dict()
    trans2ko = dict()
    for t, n_lst in t2n.items():
        for n in n_lst:
            if n in n2go:
                trans2go.setdefault(t, list())
                trans2go[t] += n2go[n]
            if n in n2k:
                path_lst = list()
                k_lst = n2k[n]
                trans2ko.setdefault(t, list())
                trans2ko[t] += k_lst
                for k in k_lst:
                    if k in k2p:
                        path_lst.extend(k2p[k])
                trans2path.setdefault(t, list())
                trans2path[t] += path_lst

    result = zip(
        ['gene2ncbi.txt', 'trans2ncbi.txt', 'gene2go.txt', 'gene2ko.txt', 'gene2path.txt',
         'trans2go.txt', 'trans2ko.txt', 'trans2path.txt'],
        [g2n, t2n, gene2go, gene2ko, gene2path,
         trans2go, trans2ko, trans2path]
    )
    for name, data in result:
        with open(name, 'w') as f:
            for k, v in data.items():
                for each_v in v:
                    f.write('{}\t{}\n'.format(k, each_v))
    return


if __name__ == '__main__':
    map_ensembl_to_go_kegg(species='hsa', tax_id='9606')
    wget.download('http://rest.kegg.jp/link/enzyme/ko', 'ko2enzyme.txt')
    from urllib.request import urlretrieve
    url = 'https://www.kegg.jp/kegg-bin/download_htext?htext=br08901&format=htext&filedir='
    urlretrieve(url, 'br08901.keg')
    wget.download('http://purl.obolibrary.org/obo/go/go-basic.obo', 'go-basic.obo')

