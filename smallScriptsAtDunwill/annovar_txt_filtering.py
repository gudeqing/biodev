import sys
import pandas as pd
files = sys.argv[1:]


for file in files:
    a = pd.read_table(file)
    f = ['pathogenic' in x.lower() for x in a['CLNSIG']]
    f0 = [x == 'reviewed_by_expert_panel' for x in a['CLNREVSTAT']]
    f1 = [x != 'intronic' for x in a['Func_refGene']]
    f2 = [x in ['BRCA1', 'BRCA2'] for x in a['Gene_refGene']]
    f3 = [any(float(y)>=0.2 for y in x.split(':')[2].split(',')) >= 0.2 for x in a['Otherinfo14']]
    f4 = [x == '.' or float(x) <= 0.01 for x in a['esp6500siv2_all']]
    f5 = [x == '.' or float(x) <= 0.01 for x in a['1000g2015aug_all']]
    f6 = [x == '.' or float(x) <= 0.01 for x in a['ExAC_ALL']]
    f7 = [x == '.' or float(x) <= 0.01 for x in a['gnomAD_exome_ALL']]
    f8 = [x == '.' or float(x) <= 0.01 for x in a['gnomAD_exome_EAS']]
    f9 = [x == '.' or float(x) <= 0.01 for x in a['ExAC_EAS']]
    m = pd.DataFrame([f1,f2,f3,f4,f5,f6,f7, f8, f9]).T
    all_match = m.apply(all, axis=1)
    p = (pd.Series(f) & pd.Series(f0) & pd.Series(f2) & pd.Series(f3)) | all_match
    r = a.loc[p]
    if r.shape[0] == 0:
        print(f'{file} is empty after filtering!')
    r.to_csv(file[:-3]+'filtered.xls', sep='\t')
