import os


def merge_picard_metrics(metric_files:list, start_marker='## METRICS', out='merged.metrics.xls', next_n=2):
    with open(out, 'w') as fw:
        find_header = False
        for each in metric_files:
            sample = os.path.basename(each).split('.', 1)[0]
            with open(each) as fr:
                while True:
                    line = fr.readline()
                    if not line:
                        break
                    if line.startswith(start_marker):
                        if not find_header:
                            header = 'sample\t' + fr.readline()
                            fw.write(header)
                            find_header = True
                        else:
                            _ = fr.readline()
                        for _ in range(next_n-1):
                            info = f'{sample}\t' + fr.readline()
                            fw.write(info)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['merge_picard_metrics'])
