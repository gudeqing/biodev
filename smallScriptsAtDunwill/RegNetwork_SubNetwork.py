def sub_network(network, query, out='sub_network', s_ind=1, t_ind=3,
                source_must_in_query=1, target_must_in_query=1, either_in_query=False):
    """
    用于转录因子分析，从母网络中提取子网络
    :param network: 母网络文件
    :param query: 基因列表文件，决定子网络的节点抽提，本程序仅提取第一列信息
    :param out: 输出的互作网络文件
    :param s_ind: 母网络中，转录因子所在列的序号
    :param g_ind: 母网络中，转录因子靶向的基因所在列的序号
    :param source_must_in_query: 是否要求子网络的source节点必须包含在query中, 0表示否，其他整数表示是
    :param target_must_in_query: 是否要求子网络的target节点必须包含在query中, 0表示否，其他整数表示是
    :param either_in_query: 指示是否只要source 或 target 包含在query中即可
    :return:
    """
    if not source_must_in_query and (not target_must_in_query):
        raise Exception('source_must_in_query and source_must_in_query cannot be both False')
    queries = set(x.strip().split()[0] for x in open(query))
    with open(network) as fr, open(out, 'w') as fw:
        for line in fr:
            lst = line.strip().split()
            source = lst[s_ind-1]
            target = lst[t_ind-1]
            judge = source in queries
            judge2 = target in queries
            if not either_in_query:
                if source_must_in_query and target_must_in_query:
                    if judge and judge2:
                        fw.write(line)
                else:
                    if source_must_in_query:
                        if judge:
                            fw.write(line)
                    if target_must_in_query:
                        if judge2:
                            fw.write(line)
            else:
                if judge or judge2:
                    fw.write(line)


if __name__ == '__main__':
    from xcmds.xcmds import xcmds
    xcmds(locals())

