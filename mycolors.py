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
