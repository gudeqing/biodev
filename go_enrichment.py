import pandas as pd
import numpy as np
from glob import glob
from statsmodels.stats.multitest import multipletests
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.godag_plot import plot_gos
from goatools.associations import read_associations
__author__ = 'gudeqing'


def enrich(gene2go:str, study:str, obo:str, population:str=None, geneid2symbol:str=None, correct='fdr_bh',
              alpha=0.05, top=20, goea_out=None, dag_out=None, dpi=300, show_gene_limit=6, only_plot_sig=False):
    """
    Go enrichment based on goatools
    :param gene2go: a file with two columns: gene_id \t go_term_id
    :param study: a file with at least one column, first column contains gene id, second columns is regulation direction
    :param obo: go-basic file download from GeneOntology
    :param population: a file with each row contains one gene; default to use all genes in gene2go file as population
    :param geneid2symbol: file with two columns: gene_id \t gene_symbol, used for DAG plot
    :param correct: pvalue adjustment method:
        Method used for testing and adjustment of pvalues. Can be either the
        full name or initial letters. Available methods are:
        - `bonferroni` : one-step correction
        - `sidak` : one-step correction
        - `holm-sidak` : step down method using Sidak adjustments
        - `holm` : step-down method using Bonferroni adjustments
        - `simes-hochberg` : step-up method  (independent)
        - `hommel` : closed method based on Simes tests (non-negative)
        - `fdr_bh` : Benjamini/Hochberg  (non-negative)
        - `fdr_by` : Benjamini/Yekutieli (negative)
        - `fdr_tsbh` : two stage fdr correction (non-negative)
        - `fdr_tsbky` : two stage fdr correction (non-negative)
    :param alpha: fdr cutoff, default 0.05
    :param top: n top go terms to plot, sorted by corrected pvalue
    :param goea_out: output enrichment result file
    :param dag_out: dag figure file
    :param dpi: resolution of image, no effect for svg
    :param show_gene_limit: the max number of gene in a node to show
    :param only_plot_sig: only plot dag for significantly enriched terms
    :return: None
    """
    if str(correct) == '3':
        correct = 'fdr_bh'
    if geneid2symbol:
        geneid2symbol = dict(x.strip().split()[:2] for x in open(geneid2symbol) if x.strip())
    else:
        geneid2symbol = dict()
    obo = GODag(obo, optional_attrs=['relationship', 'is_a'])
    gene2go = read_associations(gene2go)
    study_genes = [x.strip().split()[0] for x in open(study)]
    try:
        reg_dict = dict(x.strip().split()[:2] for x in open(study))
    except:
        reg_dict = {x.strip(): '' for x in open(study)}
    if not population:
        population = gene2go.keys()
    else:
        population = [x.strip().split()[0] for x in open(population) if x.strip()]

    goea_obj = GOEnrichmentStudy(
        population, gene2go, obo,
        propagate_counts=False, alpha=alpha,
        methods=('fdr_bh', )
    )
    keep_if = lambda r: r.ratio_in_study[0] != 0
    goea_results_all = goea_obj.run_study(study_genes, keep_if=keep_if)
    goea_out = goea_out or study + '.goea.xls'
    goea_obj.wr_tsv(goea_out, goea_results_all)

    def func(y):
        results = []
        genes = [x.strip() for x in y.split(',')]
        for gene in genes:
            tmp = [gene]
            if gene in reg_dict:
                tmp.append(reg_dict[gene])
            if gene in geneid2symbol:
                tmp.append(geneid2symbol[gene])
            results.append('|'.join(tmp))
        return ';'.join(results)
    # func = lambda y: ';'.join(x.strip()+'|'+reg_dict[x.strip()] if x.strip() in reg_dict else x.strip() for x in y.split(','))
    table = pd.read_table(goea_out, header=0, index_col=0)
    # 重新校正pvalue, 修改内容
    fdr = multipletests(table['p_uncorrected'], method=correct)[1]
    table['p_fdr_bh'] = fdr
    # 修改goea_result_all方便后续的画图
    for r, fdr in zip(goea_results_all, fdr):
        r.p_fdr_bh = fdr
    table.columns = [x if x != 'p_fdr_bh' else 'p_corrected' for x in table.columns]
    table['enrichment'] = ['e' if x <= alpha else 'p' for x in table['p_corrected']]
    table['study_items'] = table.loc[:, 'study_items'].map(func)
    # table = table.sort_values(by=['p_corrected', 'p_uncorrected'])
    table.to_csv(goea_out, header=True, index=True, sep='\t')

    # -------------------plot dag------------------------
    for each in ['BP', 'MF', 'CC']:
        if only_plot_sig:
            goea_results_sig = table[table['enrichment'] == 'e']
        else:
            goea_results_sig = table.copy()
        goea_results_sig = goea_results_sig[goea_results_sig['NS'] == each]
        if not goea_results_sig.shape[0]:
            print(f"No significant term to plot for {each} ")
            return
        if goea_results_sig.shape[0] >= top:
            goea_results_sig = goea_results_sig.iloc[:top]
        goid_subset = list(goea_results_sig.index)
        # t = obo[goid_subset[5]]
        # for k, v in t.relationship.items():
        #     print(t, k, type(v), list(v)[0].id)
        # print(dag_out[:-4]+'.'+each+dag_out[-4:])
        dag_out = dag_out or study+'.goea.dag.svg'
        plot_gos(dag_out[:-4]+'.'+each+dag_out[-4:],
                 goid_subset,  # Source GO ids, 如果分析结果里面没有包含这个节点，则他的颜色会是苍白绿色，但这里这个情况不会出现
                 obo,
                 goea_results=goea_results_all,  # use pvals for coloring:"p_{M}".format(M=goea[0].method_flds[0].fieldname)
                 # We can further configure the plot...
                 id2symbol=geneid2symbol,  # Print study gene Symbols, not GeneIDs
                 study_items=show_gene_limit,  # Only max 6 gene Symbols on GO terms
                 items_p_line=3,  # Print 3 genes per line)
                 dpi=0 if dag_out.endswith('svg') else dpi,
                 # title="Directed Graph of enriched {} terms".format(each)
                 )


def enrich_batch(study:list, gene2go=None, obo=None, geneid2symbol=None, population=None,
                 correct='fdr_bh', alpha=0.05, top=20, show_gene_limit=6, only_plot_sig=False):
    gene2go = gene2go or "/nfs2/database/Human_gene_go_kegg_annot/NewAnnot/hsa.ensembl.gene2go.txt"
    obo = obo or "/nfs2/database/Human_gene_go_kegg_annot/NewAnnot/go-basic.obo"
    geneid2symbol = geneid2symbol or "/nfs2/database/Human_gene_go_kegg_annot/NewAnnot/hsa.ensembl.id2symbol.txt"
    population = population or "/nfs2/database/Human_gene_go_kegg_annot/NewAnnot/hsa.gene.list"
    for each in study:
        enrich(gene2go=gene2go, study=each, obo=obo, population=population,
               geneid2symbol=geneid2symbol, correct=correct,
               alpha=alpha, top=top, goea_out=None, only_plot_sig=only_plot_sig,
               dag_out=None, dpi=300, show_gene_limit=show_gene_limit)


if __name__ == '__main__':
    from xcmds.xcmds import  xcmds
    xcmds(locals(), include=['enrich', 'enrich_batch'])
