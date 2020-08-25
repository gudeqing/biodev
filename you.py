from turtle import Screen, Turtle
from random import random,randint

screen = Screen()
width, height = 800, 600
screen.setup(width, height)
screen.delay(0)
screen.title("星空")
screen.bgcolor("#141852")


def draw_rose():
    t = Turtle()
    t.up()
    t.goto(100, 50)
    t.down()
    t.speed(3)
    t.penup()
    t.left(90)
    t.fd(100)
    t.pendown()
    t.right(90)
    t.fillcolor("red")
    t.begin_fill()
    t.circle(10, 180)
    t.circle(25, 110)
    t.left(50)
    t.circle(60, 45)
    t.circle(20, 170)
    t.right(24)
    t.fd(30)
    t.left(10)
    t.circle(30, 110)
    t.fd(20)
    t.left(40)
    t.circle(90, 70)
    t.circle(30, 150)
    t.right(30)
    t.fd(15)
    t.circle(80, 90)
    t.left(15)
    t.fd(45)
    t.right(165)
    t.fd(20)
    t.left(155)
    t.circle(150, 80)
    t.left(50)
    t.circle(150, 90)
    t.end_fill()
    t.left(150)
    t.circle(-90, 70)
    t.left(20)
    t.circle(75, 105)
    t.setheading(60)
    t.circle(80, 98)
    t.circle(-90, 40)
    t.left(180)
    t.circle(90, 40)
    t.circle(-80, 98)
    t.setheading(-83)
    # 画枝叶
    t.color('green')
    t.pensize(5)
    t.fd(30)
    t.left(90)
    t.fd(25)
    t.left(45)
    t.fillcolor("green")
    t.begin_fill()
    t.circle(-80, 90)
    t.right(90)
    t.circle(-80, 90)
    t.end_fill()
    t.right(135)
    t.fd(60)
    t.left(180)
    t.fd(85)
    t.left(90)
    t.fd(80)
    t.right(90)
    t.right(45)
    t.fillcolor("green")
    t.begin_fill()
    t.circle(80, 90)
    t.left(90)
    t.circle(80, 90)
    t.end_fill()
    t.left(135)
    t.fd(60)
    t.left(180)
    t.fd(60)
    t.right(90)
    t.circle(200, 60)


def draw_moon():
    t = Turtle()
    # 画月亮
    t.up()
    t.goto(-300, 250)
    t.down()
    t.fillcolor('#ffd700')
    t.begin_fill()
    t.pencolor('#ffd700')
    t.pensize(3)
    t.circle(30)
    t.end_fill()
    t.up()
    t.hideturtle()


def generate_stars():
    t = Turtle(visible=False, shape='circle')
    t.pencolor("white")
    t.fillcolor("white")
    t.penup()
    t.setheading(-90)
    t.goto(width/2,randint(-height/2,height/2))
    stars = []
    colors = ['yellow', 'white', 'pink']
    for _ in range(125):
        star = t.clone()
        s = random() /3
        star.shapesize(s,s)
        star.speed(int(s*8))
        star.setx(width/2 + randint(1,width))
        star.sety( randint(-height/2,height/2))
        star.color(colors[randint(0,2)], colors[randint(0,2)])
        star.begin_fill()
        star.showturtle()
        star.end_fill()
        stars.append(star)
    return stars


def main():
    t = Turtle()
    t.hideturtle()
    i = 0
    stars = generate_stars()
    draw_moon()
    while True:
        i += 1
        if i == 125:
            draw_rose()
        elif i == 520:
            t.up()
            t.goto(30, -190)
            t.down()
            t.color('white')
            t.write('Best wishes to 张慧', font=('STliti', 18, 'italic'), align="center")
            t.up()
        elif i == 2000:
            t.up()
            t.goto(0, 0)
            t.down()
            t.write('让眼睛休息一下@-@', font=('STliti', 20, 'italic'), align="center")
            t.hideturtle()
            t.up()
        elif i == 50000:
            t.up()
            t.goto(0, 100)
            t.down()
            t.write('你太可爱了', font=('STliti', 20, 'italic'), align="center")
            t.hideturtle()
            t.up()

        for star in stars:
            star.setx(star.xcor() - 3 * star.speed())
            star.speed(3)
            if star.xcor()<-width/2:
                star.hideturtle()
                star.setx(width/2 + randint(1,width))
                star.sety( randint(-height/2,height/2))
                star.showturtle()


if __name__ == '__main__':
    main()

# pyinstaller -F you.py --clean
