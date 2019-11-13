import os
import pandas as pd


def merge_star_alignment_stat(log_files: list, outfile=None):
    results = list()
    for logfile in log_files:
        sample = os.path.basename(logfile).split('.', 1)[0]
        with open(logfile) as fr:
            _ = [fr.readline() for i in range(5)]
            result = dict(sample=sample)
            for line in fr:
                if '|' in line:
                    desc, value = line.split('|')
                    desc = desc.strip()
                    value = value.strip()
                    result[desc] = value
            results.append(result)
    df = pd.DataFrame(results).set_index('sample')
    outfile = 'star_alignment_stat.csv' if outfile is None else outfile
    df.to_csv(outfile)
    return df


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['merge_star_alignment_stat'])
