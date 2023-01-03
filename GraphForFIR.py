# coding=utf8
import sys

import matplotlib.pyplot as plt
import asyncio
import numpy as np
import pylab
from sys import argv
import wave as wv
import scipy.io.wavfile as wv1

print(sys.argv[1])

if(sys.argv[1] == 'graph1'):
    data = []
    x = []
    print('ITS PYTHON')
    file_name = sys.argv[2]
    print(file_name)
    f = open(file_name, 'r')
    it = 0
    for i in f:
        data.append(float(i))
        x.append(it)
        it = it + (1.0 / 44100.0)

    plt.plot(x, data)
    plt.grid()
    plt.show()
elif(sys.argv[1] == 'graph2'):
    data = []
    x = []
    print('ITS PYTHON')
    file_name = sys.argv[2]
    print(file_name)
    f = open(file_name, 'r')
    it = 0
    for i in f:
        data.append(float(i))
        x.append(it)
        it = it + (1.0/44100)

    plt.plot(x, data)
    plt.grid()
    plt.show()
elif (sys.argv[1] == 'graph3'):
    data = []
    data1 = []
    data2 = []
    x = []
    print('ITS PYTHON')
    file_name = sys.argv[2]
    print(file_name)
    f = open(file_name, 'r')
    it = 0.0
    for i in f:
        data.append(float(i))
        x.append(it)
        it = it + (1.0/44100.0)
    data1 = data[:int(len(data)/2)]
    data2 = data[int(len(data)/2):]
    x = x[:int(len(data)/2)]
    plt.plot(x, data1)
    plt.plot(x, data2)
    plt.grid()
    plt.show()
elif (sys.argv[1] == 'graph4'):
    data1 = []
    data2 = []
    x1 = []
    x2 = []
    print('ITS PYTHON')
    file_name = sys.argv[2]
    lenData = int(sys.argv[3]) - 4
    print(file_name)
    f = open(file_name)
    it = 0
    it1 = 0
    it2 = 0
    for line in f:
        if(it <= lenData/2):
            data1.append(float(line))
            it1 += (1 / 44100)
            x1.append(it1)
        else:
            data2.append(float(line))
            it2 += (1 / 44100)
            x2.append(it2)
        it = it + 1
    plt.plot(x1, data1)
    plt.plot(x2, data2)
    plt.grid()
    plt.show()
elif(sys.argv[1] == 'readWav'):
    file_wave = wv.open('C:\\Users\\Пользователь\\source\\repos\\RpojectOOOSTC1\\RpojectOOOSTC1\\file_read.wav', 'rb')
    print("file has been uploading")
    mass_signal = list(file_wave.readframes(file_wave.getnframes()))
    max_mass_signal = max(mass_signal)

    for i in range(0, len(mass_signal)):
        mass_signal[i] = mass_signal[i]

    file_read_wav = open('C:\\Users\\Пользователь\\source\\repos\\RpojectOOOSTC1\\RpojectOOOSTC1\\file_read.txt', 'w')
    file_read_wav.write(str(file_wave.getframerate()))
    file_read_wav.write('\n')
    file_read_wav.write(str(file_wave.getnframes()))
    file_read_wav.write('\n')

    for i in range(0, len(mass_signal)):
        if(i % 4 == 0):
            file_read_wav.write(str(mass_signal[i]) + '\n')

elif(sys.argv[1] == 'writeWav'):
    data = []
    file_name = str(sys.argv[2])
    print(file_name + '.txt')
    with open(file_name + '.txt', 'r') as f:
        for line in f:
            data.append([int(x) for x in line.split()])
    wv1.write(file_name + '.wav', 44100, np.int16(data))