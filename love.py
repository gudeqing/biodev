from turtle import *
from random import random,randint
import time

screen = Screen()
width ,height = 800,600
screen.setup(width,height)
screen.delay(0)
screen.title("星空")
screen.bgcolor("black")


def draw_rose():
    import turtle
    turtle.goto(100, 50)
    turtle.speed(1)
    turtle.penup()
    turtle.left(90)
    turtle.fd(100)
    turtle.pendown()
    turtle.right(90)
    turtle.fillcolor("red")
    turtle.begin_fill()
    turtle.circle(10, 180)
    turtle.circle(25, 110)
    turtle.left(50)
    turtle.circle(60, 45)
    turtle.circle(20, 170)
    turtle.right(24)
    turtle.fd(30)
    turtle.left(10)
    turtle.circle(30, 110)
    turtle.fd(20)
    turtle.left(40)
    turtle.circle(90, 70)
    turtle.circle(30, 150)
    turtle.right(30)
    turtle.fd(15)
    turtle.circle(80, 90)
    turtle.left(15)
    turtle.fd(45)
    turtle.right(165)
    turtle.fd(20)
    turtle.left(155)
    turtle.circle(150, 80)
    turtle.left(50)
    turtle.circle(150, 90)
    turtle.end_fill()
    turtle.left(150)
    turtle.circle(-90, 70)
    turtle.left(20)
    turtle.circle(75, 105)
    turtle.setheading(60)
    turtle.circle(80, 98)
    turtle.circle(-90, 40)
    turtle.left(180)
    turtle.circle(90, 40)
    turtle.circle(-80, 98)
    turtle.setheading(-83)
    # 画枝叶
    turtle.color('green')
    turtle.pensize(3)
    turtle.fd(30)
    turtle.left(90)
    turtle.fd(25)
    turtle.left(45)
    turtle.fillcolor("green")
    turtle.begin_fill()
    turtle.circle(-80, 90)
    turtle.right(90)
    turtle.circle(-80, 90)
    turtle.end_fill()
    turtle.right(135)
    turtle.fd(60)
    turtle.left(180)
    turtle.fd(85)
    turtle.left(90)
    turtle.fd(80)
    turtle.right(90)
    turtle.right(45)
    turtle.fillcolor("green")
    turtle.begin_fill()
    turtle.circle(80, 90)
    turtle.left(90)
    turtle.circle(80, 90)
    turtle.end_fill()
    turtle.left(135)
    turtle.fd(60)
    turtle.left(180)
    turtle.fd(60)
    turtle.right(90)
    turtle.circle(200, 60)


t = Turtle(visible = False,shape='circle')
t.pencolor("white")
t.fillcolor("white")
t.penup()
t.setheading(-90)
t.goto(width/2,randint(-height/2,height/2))
stars = []
colors = ['yellow', 'white', 'pink']

# 画月亮
t.goto(-300, 270)
t.fillcolor('#ffd700')#填充色
t.begin_fill()
t.pencolor('#ffd700')
t.pensize(3)
t.circle(30)
t.end_fill()

for i in range(125):
    star = t.clone()
    s = random() /3
    star.shapesize(s,s)
    star.speed(int(s*10))
    star.setx(width/2 + randint(1,width))
    star.sety( randint(-height/2,height/2))
    star.color(colors[randint(0,2)], colors[randint(0,2)])
    star.begin_fill()
    star.showturtle()
    star.end_fill()
    stars.append(star)

# draw_rose()
i = 0
while True:
    i += 1
    if i == 125:
        draw_rose()
    elif i == 520:
        star.goto(30, -190)
        star.down()
        star.write('Best wishes to 张慧', font=('STliti', 18, 'italic'), align="center")
        star.up()
    elif i == 2000:
        star.goto(0, 0)
        star.write('让眼睛休息一下@-@', font=('STliti', 18, 'italic'), align="center")

    for star in stars:
        star.setx(star.xcor() - 3 * star.speed())
        if star.xcor()<-width/2:
            star.hideturtle()
            star.setx(width/2 + randint(1,width))
            star.sety( randint(-height/2,height/2))
            star.showturtle()

