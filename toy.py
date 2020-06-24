import turtle
from turtle import Screen
import random


def draw_heart(t, pos, size=100):
    """
    设想一个正菱形，分别以其中两条相邻的边作为直径画半圆，最后由两条边和两个半圆构成心形图案
    """
    t.up()
    t.goto(pos)
    t.down()
    # t.showturtle()

    # 设置线条大小
    t.pensize(2)
    # 设置颜色，线条颜色为red，填充色为粉色
    t.color("red", 'pink')
    t.begin_fill()
    # 画边 逆时针转45° 沿着当前方向走200
    t.left(45)
    t.fd(size)

    # 画半圆
    """
    circle(radius, extent)的参数说明
    radius ：弧形半径
        当radius值为正数时，圆心在当前位置 / 小海龟左侧。
        当radius值为负数时，圆心在当前位置 / 小海龟右侧。
    extent ：弧形角度。当无该参数或参数为None时，绘制整个圆形
        当extent值为正数时，顺小海龟当前方向绘制。
        当extent值为负数时，逆小海龟当前方向绘制。
    """
    t.circle(size/2, 180)
    # 顺时针转向90°
    t.right(90)
    # 画半圆
    t.circle(size/2, 180)
    # 画边
    t.fd(size)
    # 填充颜色
    t.end_fill()


def draw_tree(n, branch_ratio=0.77, angle=30):
    pos_lst = []
    heart_size = []
    # 下面是一个递归函数画树
    def tree(size):
        size = size * random.randrange(80, 100) * 0.01
        SIZE_TREE = 15 * random.randrange(5, 15)*0.1
        turtle.pensize(size/10)
        if size > SIZE_TREE:
            turtle.forward(size)
            turtle.right(angle)
            tree(size * branch_ratio)
            # 左边
            turtle.left(angle*2)
            tree(size * branch_ratio)

            # 回到之前的树枝
            turtle.right(angle)
            # 给最后的树枝画绿色
            if size / 2 <= SIZE_TREE:
                turtle.color('green')
            else:
                turtle.color('brown')
            turtle.backward(size)
        else:
            pos = turtle.pos()
            if pos not in pos_lst:
                pos_lst.append(pos)
                heart_size.append(size)
    tree(n)
    return pos_lst, heart_size


def main():
    # screensize 参数分别为画布的宽(单位像素),高,背景颜色
    canvas = Screen()
    canvas.screensize(800, 600, bg='wheat')

    # 画树
    turtle.penup()
    turtle.goto(0,-120)
    turtle.speed(0)
    turtle.hideturtle()
    turtle.left(75)
    turtle.backward(100)
    # turtle.showturtle()
    turtle.pendown()
    turtle.color('brown')
    pos_lst, heart_size = draw_tree(150)

    # 画心
    for pos, size in zip(pos_lst, heart_size):
        draw_heart(turtle, pos, size=size*0.8)

    # 落款
    turtle.hideturtle()
    turtle.color('black','pink')
    turtle.up()
    turtle.goto(120,-160)
    turtle.speed(5)
    turtle.down()
    turtle.write('Best wishes to ？',font=('STliti',15, 'italic'), align="center")
    turtle.up()
    turtle.goto(120,-190)
    turtle.down()
    turtle.write('2020.05.20 from gdq', font=('STliti',12, 'italic'), align="center")

    # 点击图片退出
    canvas.exitonclick()


if __name__ == '__main__':
    main()