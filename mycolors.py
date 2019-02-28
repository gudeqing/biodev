import colorlover as cl
from IPython.display import HTML
from colorsys import hls_to_rgb
import numpy as np


def get_color_pool(n):
    if n <= 12:
        return cl.scales['12']['qual']['Paired']
    color_pool = []
    for i in np.arange(60., 360., 360. / n):
        hue = i / 300.
        rand_num = np.random.random_sample()
        lightness = (50 + rand_num * 10) / 100.
        saturation = (90 + rand_num * 10) / 100.
        rgb = hls_to_rgb(hue, lightness, saturation)
        color_pool.append(tuple([int(x * 255) for x in rgb]))
    return cl.to_rgb(color_pool)


# 随机生成100不同的颜色
n = 100
p = [tuple(x) for x in np.random.choice(range(10, 250), (n,3))]
p = list(set(p))
p.sort(key=lambda x:x)
color_pool = cl.to_rgb(p[::-1])

# 获取12个漂亮的颜色
color_pool = cl.scales['12']['qual']['Paired']





def ucsc2gencode(infile):
    return dict(x.strip().split()[:2] for x in open(infile))


def replace_ucsc_with_gencode(infile, outfile, convert_file='GRCh38_UCSC2gencode.txt'):
    convert_dict = ucsc2gencode(convert_file)
    with open(infile) as fr, open(outfile, 'w') as fw:
        for line in fr:
            chr_name = line.split()[0]
            if chr_name in convert_dict:
                line = line.replace(chr_name, convert_dict[chr_name], 1)
                fw.write(line)
            else:
                fw.write(line)

def filter_bed_of_target_chr(infile, outfile, chr_list_file='primary_assembly.chr.list'):
    chr_list = [x.strip() for x in  open(chr_list_file)]
    with open(infile) as fr, open(outfile, 'w') as fw:
        for line in fr:
            chr_name = line.split()[0]
            if chr_name in chr_list:
                fw.write(line)
            else:
                pass
def get_target_gene_matrix(infile, outfile, target, filter_out_target=False):
    with open(infile) as fr, open(outfile, 'w') as fw, open(target) as fr2:
        target_list = set(x.strip().split('.')[0] for x in fr2)
        fw.write(fr.readline())
        for line in fr:
            gene = line.split('\t', 1)[0].split('_')[0].split('.', 1)[0]
            if filter_out_target:
                if gene not in target_list:
                    fw.write(line)
            else:
                if gene in target_list:
                    fw.write(line)

if __name__ == '__main__':
    def introduce_command(func):
        import argparse
        import inspect
        import json
        import time
        if isinstance(func, type):
            description = func.__init__.__doc__
        else:
            description = func.__doc__
        parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
        func_args = inspect.getfullargspec(func)
        arg_names = func_args.args
        arg_defaults = func_args.defaults
        arg_defaults = ['None']*(len(arg_names) - len(arg_defaults)) + list(arg_defaults)
        for arg, value in zip(arg_names, arg_defaults):
            if arg == 'self':
                continue
            if value == 'None':
                parser.add_argument('-'+arg, required=True, metavar=arg)
            elif type(value) == bool:
                if value:
                    parser.add_argument('--'+arg, action="store_false", help='default: True')
                else:
                    parser.add_argument('--'+arg, action="store_true", help='default: False')
            elif value is None:
                parser.add_argument('-' + arg, default=value, metavar='Default:' + str(value), )
            else:
                parser.add_argument('-' + arg, default=value, type=type(value), metavar='Default:' + str(value), )
        if func_args.varargs is not None:
            print("warning: *varargs is not supported, and will be neglected! ")
        if func_args.varkw is not None:
            print("warning: **keywords args is not supported, and will be neglected! ")
        args = parser.parse_args().__dict__
        with open("Argument_detail.json", 'w') as f:
            json.dump(args, f, indent=2, sort_keys=True)
        start = time.time()
        func(**args)
        print("total time: {}s".format(time.time() - start))

    introduce_command(get_target_gene_matrix)







