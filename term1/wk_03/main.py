import matplotlib.pyplot
import matplotlib.pyplot as plt
import numpy as np

f = open('/home/starman/CLionProjects/EulersMethod/data.txt', 'r')
conv = {
    0: lambda x: float(x),
    1: lambda x: float(x),
    2: lambda x: float(x),
    3: lambda x: float(x),
}
arr = np.loadtxt(f, delimiter=" ", converters=conv)

t, x, energy, v = [], [], [], []
for i in arr:
    t.append(i[0])
    x.append(i[1])
    energy.append(i[2])
    v.append(i[3])

lines = 1
columns = 3
fig, ax = plt.subplots(lines, columns)
plt.subplot(lines, columns, 1)
plt.plot(t, x)
plt.title('dt = 0.01')
plt.subplot(lines, columns, 2)
plt.plot(t, energy)
plt.title('energy')
plt.subplot(lines, columns, 3)
plt.plot(x, v)
plt.title('x-v')

plt.show()
