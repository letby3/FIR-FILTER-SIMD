import matplotlib.pyplot as plt
import asyncio


a = 128
data = []
x = []
xk = []
y1 = []
y2 = []
y3 = []
with open("C:\\Users\\Пользователь\\source\\repos\\RpojectOOOSTC1\\RpojectOOOSTC1\\graph1.txt", 'r') as f:
    for line in f:
        data.append([float(x) for x in line.split()])

for i in range(0, a):
    xk.insert(i, i)
    x.insert(i, data[i])
    #print(i)
for i in range(a, a*2):
    y1.insert(i, data[i])
    #print(i)
for i in range(0, 18):
    it = 0
    for j in range(a*i + a*2, a*i + a*3):
        y3.insert(it, data[j])
        it+=1
    y2.insert(i, y3)
    y3 = []
'''
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
'''
plt.plot(x, y1)
plt.grid()
#plt.figure(1)
plt.plot(x, y2[2])
plt.plot(x, y2[3])
plt.show()
