#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <conio.h>
#include <vector>
#include <complex>
#include <time.h>
#include "FIR_FILTER.h"

using namespace std;

ofstream out;

typedef complex<double> base;

void FIR_Filter(long float x[], long float y[], const int size_n);

void fft(vector<base>& a, bool invert);


int main(){
    srand(time(NULL));
    out.open("graph.txt"); //Файл вывода 
    long float N = 128; //Точек отчеста входного сигнала 
    long float fs = 10; //Частота дискретизации (в сек)
    long float X_MAX = N / fs - 1/fs; //Длина сигнала в сек
    long float *x = new long float[1001]; //Входной сигнал
    long float *y = new long float[1001]; //Выходной сигнал
    long float* yY = new long float[1001];
    int it = 0; //верный друг, верный счётчик
    vector<base> a, b; // ВХОД и ВЫХОД Фурье

    for (long float i = 0, j = 0; i <= X_MAX;i += 1.0/fs, j++) {
        out << i << endl; //Вывод координаты ВХ сигнала в сек
    }

    for (long float i = 0; i <= X_MAX; i += 1.0 / fs, it++) {
        long float X = ((long float)(abs(rand() % 12)) / 100.0 - 0.01) * cos(2 * M_PI * 0.5 * i); // Rand-ый sin-ый входной сигнал
        out << X << endl;
        x[it] = X;
        a.push_back(X);
    }

    fft(a, 0); //Фурье входного сигнала (Complex)

    FIR_Filter(x, y, N); // КИХ фильтр нижних частот

    for (int i = 0; i < N;i++) {
        out << y[i] << endl; // Вывод выходного сигнала 
        b.push_back(y[i]);
    }

    fft(b, 0); //Фурье выходного сигнала (Complex)

    for (complex<double> i : a) {
        out << i.real() << endl; //Вывод входного Фурье
    }

    for (complex<double> i : b) {
        out << i.real() << endl; //Вывод выходного Фурье
    }

    system("PAUSE");
    return 0;
}
