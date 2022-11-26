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

long float* FIR_Filter(long float x[], long float y[], const int size_n);

void fft(vector<base>& a, bool invert);


int main(){
    srand(time(NULL));
    out.open("graph.txt");
    long float N = 128;
    long float fs = 10;
    long float X_MAX = N / fs - 1/fs;
    long float *x = new long float[1001];
    long float *y = new long float[1001];
    int it = 0;
    vector<base> a, b;

    for (long float i = 0, j = 0; i <= X_MAX;i += 1.0/fs, j++) {
        out << i << endl;
    }
    for (long float i = 0; i <= X_MAX; i += 1.0 / fs, it++) {
        long float X = ((long float)(abs(rand() % 12)) / 100.0 - 0.01) * cos(2 * M_PI * 0.5 * i);
        out << X << endl;
        x[it] = X;
        a.push_back(X);
    }

    fft(a, 0);
    y = FIR_Filter(x, y, N);
    
    for (int i = 0; i < N;i++) {
        out << y[i] << endl;
        b.push_back(y[i]);
    }

    fft(b, 0);

    for (complex<double> i : a) {
        out << i.real() << endl;
    }

    for (complex<double> i : b) {
        out << i.real() << endl;
    }

    system("PAUSE");
    return 0;
}
