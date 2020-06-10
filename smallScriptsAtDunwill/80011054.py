
def parse(txt, out='out.xls'):
    with open(txt) as f, open(out, 'w') as f2:
        for line in f:
            if not line.startswith('>'):
                lst = line.split()
                ind = [i for i, x in enumerate(lst) if x.endswith('_at')][0]
                probe = lst[ind]
                gene_name = lst[ind-1]
                cell = ' '.join(lst[:ind-1])
                entrez = lst[ind+1]
                f2.write(f'{cell}\t{gene_name}\t{probe}\t{entrez}\n')
            else:
                f2.write(line)

if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['parse'])



