import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime as dt
from collections import deque
from matplotlib.dates import date2num
from matplotlib import animation, rc
from matplotlib.animation import PillowWriter



#plt.style.use('_mpl-gallery')

# make data
x = np.sin(np.linspace(0, 10, 1000))
y = np.sin(3*np.linspace(0, 10, 1000)+np.pi/2)

# plot
lines = 2
columns = 2
fig, axs = plt.subplots(lines, columns)


#ax.plot(x, y, linewidth=2.0)

for i in range(0, lines):
    for j in range(0, columns):
        axs[i, j].set(xlim=(-3, 3), xticks=np.arange(-3, 3),
            ylim=(-3, 3), yticks=np.arange(-3, 3))

    axs[1, 1].set(xlim=(-1, 1), xticks=np.arange(-3, 3),
            ylim=(-3, 3), yticks=np.arange(-3, 3))

dots = []
hisDots = []
hisX = []
hisY = []
def foo(x, y):
    for i in x:
        for j in y:
            dots.append([x, y])
    for i in dots:
        counter = -1
        for j in dots:
            if (i == j):
                counter+=1
        hisX.append(i[0])
        hisY.append(counter)
    return 0
axs[1, 0].bar(x, y)
axs[1, 1].hist(y, bins = 10000)

# animation
line1, = axs[0, 0].plot([], [], lw=3)
line2, = axs[0, 1].plot([], [], lw=3)
line3, = axs[1, 0].plot([], [], lw=3)
line = [line1, line2, line3]

def init():
    for lineq in line:
        lineq.set_data([], [])
    return line
def animate(i):

    x1 = np.sin(np.linspace(0, 10, 1000))
    y1 = np.sin(3*np.linspace(0, 10, 1000) + 2 * np.pi * (0.01 * i))
    line[0].set_data(x1, y1)

    x2 = np.sin(np.linspace(0, 10, 1000))
    y2 = np.sin(10 * np.linspace(0, 10, 1000) + 2 * np.pi * (0.01 * i))
    line[1].set_data(x2, y2)

    return line

anim = animation.FuncAnimation(fig, animate, init_func=init,
               frames=200, interval=20, blit=True)

# save gif
#rc('animation', html='html5')
#anim.save('anime.gif', writer='pillow', fps=60)

plt.show()

