#! /data/users/dqgu/anaconda3/bin/python


def run(files:list, exp=None, out=None, has_header=False,
        intersect_only=True, intersect_xoy=1, union_only=False):
    """
    根据文件内容构建集合，并按指定规则进行运算，默认计算所有集合的交集
    :param files: 当仅提供一个文件时，文件的各列被当作是集合，集合的元素是单元格的内容；
    提供多个文件时，对每个文件的内容被一个集合，集合的元素为一整行。
    :param exp: 表达式，字符串的形式，如's1-s2'表示第一个集合减去第二个集合，集合顺序与文件提供的顺序一一对应
    :param out: 指定集合运算结果的文件名
    :param has_header: 指定文件是否包含header，默认无，如有header，header不参与计算
    :param intersect_only: 如提供，不考虑exp指定的运算，而是计算所有集合的交集，即交集结果的所有元素在集合中出现的频数等于集合数。
    :param intersect_xoy: 如提供，不考虑exp指定的运算，而是计算所有集合的交集，而且输出交集结果的元素
    在所有集合中出现的频数大于或等于该参数指定的阈值。
    :param union_only: 计算各个集合的并集
    :return: None
    """
    set_number = len(files)
    if len(files) >= 2:
        for ind, each in enumerate(files, start=1):
            exec('s{}=set(open("{}").readlines())'.format(ind, each))
    else:
        import pandas as pd
        table = pd.read_table(files[0], header=0 if has_header else None)
        set_number = table.shape[1]
        for i in range(table.shape[1]):
            exec('s{}=set(table.iloc[:, {}])'.format(i+1, i))

    result = list()
    if exp:
        print("do as you say in exp")
        result = eval(exp)
    elif intersect_xoy > 1:
        print('do intersect_xoy')
        union = eval('|'.join(['s'+str(x) for x in range(1, set_number+1)]))
        result = set()
        for each in union:
            varspace = dict(locals())
            in_times = sum(eval("each in s{}".format(x), varspace) for x in range(1, set_number+1))
            if in_times >= intersect_xoy:
                result.add(each)
    elif union_only:
        print('do union only')
        result = eval('|'.join(['s'+str(x) for x in range(1, set_number+1)]))
    elif intersect_only:
        print('do intersect only')
        result = eval('&'.join(['s'+str(x) for x in range(1, set_number+1)]))
    if not result:
        print('result is empty!')
    else:
        print('result size: {}'.format(len(result)))
    with open(out or 'result.list', 'w') as f:
        _ = [f.write(x) for x in result]


if __name__ == '__main__':
    from xcmds.xcmds import xcmds
    xcmds(locals(), include=['run'])
