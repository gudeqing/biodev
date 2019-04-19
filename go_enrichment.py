import pandas as pd
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.godag_plot import plot_gos
from goatools.associations import read_associations
__author__ = 'gudeqing'


def enrich(gene2go:str, study:str, obo:str, population:str=None, geneid2symbol:str=None, correct=('fdr_bh', ),
              alpha=0.05, top=20, goea_out='goea.xls', dag_out='goea.svg', dpi=300, show_gene_limit=6):
    """
    Go enrichment based on goatools
    :param gene2go: a file with two columns: gene_id \t go_term_id
    :param study: a file with at least one column, first column contains gene id, second columns is regulation direction
    :param obo: go-basic file download from GeneOntology
    :param population: a file with each row contains one gene; default to use all genes in gene2go file as population
    :param geneid2symbol: file with two columns: gene_id \t gene_symbol, used for DAG plot
    :param correct: pvalue adjustment method
    :param alpha: fdr cutoff, default 0.05
    :param top: n top go terms to plot, sorted by corrected pvalue
    :param goea_out: output enrichment result file
    :param dag_out: dag figure file
    :param dpi: resolution of image, no effect for svg
    :param show_gene_limit: the max number of gene in a node to show
    :return: None
    """
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
        reg_dict = dict()
    if not population:
        population = gene2go.keys()
    else:
        population = [x.strip().split()[0] for x in open(population) if x.strip()]

    goea_obj = GOEnrichmentStudy(
        population, gene2go, obo,
        propagate_counts=False, alpha=alpha,
        methods=correct
    )
    keep_if = lambda r: r.ratio_in_study[0] != 0
    goea_results_all = goea_obj.run_study(study_genes, keep_if=keep_if)
    goea_obj.wr_tsv(goea_out, goea_results_all)
    if reg_dict:
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
        table['study_items'] = table.loc[:, 'study_items'].map(func)
        table.to_csv(goea_out, header=True, index=True, sep='\t')

    # -------------------plot dag------------------------
    for each in ['BP', 'MF', 'CC']:
        goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < alpha and r.NS == each]
        if not goea_results_sig:
            print("No significant term to plot, and exit now")
            return
        if len(goea_results_sig) >= top:
            goea_results_sig = goea_results_sig[:top]
        goid_subset = [r.GO for r in goea_results_sig]
        # t = obo[goid_subset[5]]
        # for k, v in t.relationship.items():
        #     print(t, k, type(v), list(v)[0].id)
        # print(dag_out[:-4]+'.'+each+dag_out[-4:])
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


if __name__ == '__main__':
    class Func2Command(object):
        def __init__(self, callable_dict):
            self.call(callable_dict)

        @staticmethod
        def introduce_command(func):
            import argparse
            import inspect
            import json
            import time
            import sys
            if isinstance(func, type):
                description = func.__init__.__doc__
            else:
                description = func.__doc__
            if '-h' not in sys.argv or '--help' in sys.argv or '-help' in sys.argv:
                description = None
            if description:
                _ = [print(x.strip()) for x in description.split('\n') if x.strip()]
                parser = argparse.ArgumentParser(add_help=False)
            else:
                parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=description)
            func_args = inspect.getfullargspec(func)
            arg_names = func_args.args
            if not arg_names:
                func()
                return
            arg_defaults = func_args.defaults
            if not arg_defaults:
                arg_defaults = []
            arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
            sig = inspect.signature(func)
            for arg, value in zip(arg_names, arg_defaults):
                arg_type = sig.parameters[arg].annotation
                if arg == 'self':
                    continue
                if value == 'None':
                    if arg_type in [list, tuple, set]:
                        parser.add_argument('-' + arg, nargs='+', required=True, metavar=arg)
                    elif arg_type in [int, str, float]:
                        parser.add_argument('-' + arg, type=arg_type, required=True, metavar=arg)
                    else:
                        parser.add_argument('-'+arg, required=True, metavar=arg)
                elif type(value) == bool:
                    if value:
                        parser.add_argument('--'+arg, action="store_false", help='default: True')
                    else:
                        parser.add_argument('--'+arg, action="store_true", help='default: False')
                elif value is None:
                    if arg_type in [list, tuple, set]:
                        parser.add_argument('-' + arg, nargs='+', required=False, metavar=arg)
                    elif arg_type in [int, str, float]:
                        parser.add_argument('-' + arg, type=arg_type, required=False, metavar=arg)
                    else:
                        parser.add_argument('-'+arg, required=False, metavar='Default:' + str(value))
                else:
                    if arg_type in [list, tuple, set] or (type(value) in [list, tuple, set]):
                        default_value = ' '.join(str(x) for x in value)
                        if type(value) in [list, tuple]:
                            one_value = value[0]
                        else:
                            one_value = value.pop()
                        parser.add_argument('-' + arg, default=value, nargs='+', type=type(one_value),
                                            metavar='Default:'+default_value, )
                    else:
                        parser.add_argument('-' + arg, default=value, type=type(value),
                                            metavar='Default:' + str(value), )
            if func_args.varargs is not None:
                print("warning: *varargs is not supported, and will be neglected! ")
            if func_args.varkw is not None:
                print("warning: **keywords args is not supported, and will be neglected! ")
            args = parser.parse_args().__dict__
            try:
                with open("Argument_detail.json", 'w') as f:
                    json.dump(args, f, indent=2, sort_keys=True)
            except IOError:
                print('Current Directory is not writable, thus argument log is not written !')
            start = time.time()
            func(**args)
            print("total time: {}s".format(time.time() - start))

        def call(self, callable_dict):
            import sys
            excludes = ['introduce_command', 'Func2Command']
            _ = [callable_dict.pop(x) for x in excludes if x in callable_dict]
            if len(callable_dict) >= 2:
                if len(sys.argv) <= 1:
                    print("The tool has the following sub-commands: ")
                    _ = [print(x) for x in callable_dict]
                    return
                sub_cmd = sys.argv[1]
                sys.argv.remove(sub_cmd)

                if sub_cmd in callable_dict:
                    self.introduce_command(callable_dict[sub_cmd])
                else:
                    print('sub-command: {} is not defined'.format(sub_cmd))
            elif len(callable_dict) == 1:
                self.introduce_command(callable_dict.pop(list(callable_dict.keys())[0]))

    callable_dict = {x: y for x, y in locals().items() if callable(y)}
    _ = [callable_dict.pop(x) for x in {'Func2Command'} if x in callable_dict]
    Func2Command(callable_dict)
