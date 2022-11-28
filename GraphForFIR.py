import matplotlib.pyplot as plt
import asyncio
import numpy as np
import pylab

a = 128
data = []
data_stat = []
x_stat = []
y_stat = []
y_stat_simd = []
x = []
xk = []
y1 = []
y2 = []
y3 = []
with open("C:\\Users\\letby\\source\\repos\\VS studio Project 1\\graph1.txt", 'r') as f:
    for line in f:
        data.append([float(x) for x in line.split()])
with open("C:\\Users\\letby\\source\\repos\\VS studio Project 1\\stat_out.txt", 'r') as f_stat:
    for line in f_stat:
        data_stat.append([float(x) for x in line.split()])

it = 0

for i in range(0, len(data_stat)):
    y_stat.insert(it, data_stat[i][0])
    y_stat_simd.insert(it, data_stat[i][1])
    x_stat.insert(it, data_stat[i][2])
    it+=1

print(y_stat)

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
plt.plot(x, y1)
plt.grid()
#plt.figure(1)
plt.plot(x, y2[2])
plt.plot(x, y2[3])
plt.show()

for i in range(0, 2, 18):
    plt.figure(i / 2)
    fig, ax = plt.subplots()
    ax.plot(x, y1)
    ax.plot(x, y2[i])
    ax.plot(x, y2[i + 1])
    ax.grid()
    ax.set_xlabel('Time')
    ax.set_ylabel('Amplitude')
    plt.show()
'''
pylab.subplot (2, 2, 1)
pylab.plot(x, y1)
pylab.plot(x, y2[0])
pylab.plot(x, y2[1])
pylab.ylabel('amplitude')
pylab.xlabel('time')
pylab.subplot (2, 2, 2)
pylab.plot(x, y1)
pylab.plot(x, y2[2])
pylab.plot(x, y2[3])
pylab.ylabel('amplitude')
pylab.xlabel('time')
pylab.subplot (2, 2, 3)
pylab.plot(x, y1)
pylab.plot(x, y2[4])
pylab.plot(x, y2[5])
pylab.ylabel('amplitude')
pylab.xlabel('time')
pylab.subplot (2, 2, 4)
pylab.plot(x, y1)
pylab.plot(x, y2[6])
pylab.plot(x, y2[7])
pylab.ylabel('amplitude')
pylab.xlabel('time')
pylab.show()
pylab.subplot (2, 2, 1)
pylab.plot(x, y1)
pylab.plot(x, y2[8])
pylab.plot(x, y2[9])
pylab.ylabel('amplitude')
pylab.xlabel('time')
pylab.subplot (2, 2, 2)
pylab.plot(x, y1)
pylab.plot(x, y2[10])
pylab.plot(x, y2[11])
pylab.ylabel('amplitude')
pylab.xlabel('time')
pylab.subplot (2, 2, 3)
pylab.plot(x, y1)
pylab.plot(x, y2[12])
pylab.plot(x, y2[13])
pylab.ylabel('amplitude')
pylab.xlabel('time')
pylab.subplot (2, 2, 4)
pylab.plot(x, y1)
pylab.plot(x, y2[14])
pylab.plot(x, y2[15])
pylab.ylabel('amplitude')
pylab.xlabel('time')
pylab.show()
pylab.subplot (1, 1, 1)
pylab.plot(x, y1)
pylab.plot(x, y2[16])
pylab.plot(x, y2[17])
pylab.ylabel('amplitude')
pylab.xlabel('time')
pylab.show()
pylab.subplot (1, 1, 1)
pylab.plot(x_stat, y_stat)
pylab.plot(x_stat, y_stat_simd)
pylab.ylabel('time, mks')
pylab.xlabel('len_filter)')
pylab.show()