import matplotlib.pyplot as plt
import asyncio


a = 128
data = []
x = []
xk = []
y1 = []
y2 = []
y3 = []
y4 = []
with open("C:\\Users\\letby\\source\\repos\\VS studio Project 1\\VS studio Project 1\\graph.txt", 'r') as f:
    for line in f:
        data.append([float(x) for x in line.split()])

for i in range(0, a):
    xk.insert(i, i)
    x.insert(i, data[i])
for i in range(a, 2 * a):
    y1.insert(i, data[i])
for i in range(2*a, 3 * a):
    y2.insert(i, data[i])
for i in range(3 * a, 4 * a):
    y3.insert(i, data[i])
for i in range(4 * a, 5 * a):
    y4.insert(i, data[i])

plt.figure(0)
plt.plot(x, y1)
plt.grid()
#plt.figure(1)
plt.plot(x, y2)
plt.grid()
plt.figure(1)
plt.plot(xk, y3)
plt.grid()
plt.plot(xk, y4)
plt.grid()
plt.show()


