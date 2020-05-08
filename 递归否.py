def fact_iter(num,product=1):
    if num==1:
        yield product
    yield fact_iter(num-1,num*product)


import types


def tramp(gen, *args, **kwargs):
    g = gen(*args, **kwargs)
    while isinstance(g, types.GeneratorType):
        g = g.__next__()
    return g

print(tramp(fact_iter, 100000))





def max_distance_to_root(tree, leaf):
    path = []

    def find_root(tree, leaf):
        """"
        tree = dict(a='b', b=['c', 'f'], c=['d', 'g'], d=['e', 'r'])
                  f   e
                 /   /
        a - b - c - d
                 \   \
                  g   r
        """
        nonlocal path
        root = leaf
        if root in tree:
            for root in tree[root]:
                path.append(root)
                yield from find_root(tree, root)
        else:
            # 除了已经返回过的节点，每次返回是最深的节点，即距离当前leaf最远的节点
            yield root, len(path), path
    return next(find_root(tree, leaf))
