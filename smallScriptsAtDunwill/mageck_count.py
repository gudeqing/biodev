import pandas as pd


def generate_cmd(fastq_info, library, control=None, norm='median', out_prefix='Result', discard_read2=True):
    df = pd.read_csv(fastq_info, header=None, index_col=0, sep=None, engine='python')
    # mageck count
    cmd = 'mageck count '
    cmd += '-l {} '.format(library)
    cmd += '--norm-method {} '.format(norm)
    cmd += '--sample-label {} '.format(','.join(df.index))
    cmd += '-n {} '.format(out_prefix)
    if control is not None:
        cmd += '--control-sgrna {} '.format(control)

    cmd += '--fastq {} '.format(' '.join(df[1]).replace(';', ','))
    if df.shape[1] >=2:
        if not discard_read2:
            cmd += '--fastq-2 {} '.format(' '.join(df[2]).replace(';', ','))
            cmd += '--count-pair True '

    with open('mageck.count.sh', 'w') as f:
        f.write(cmd + '\n')
    return cmd


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), exclude=['pd'])
