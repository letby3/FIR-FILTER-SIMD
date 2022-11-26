#include <iostream>
#include <cmath>
#include <fstream>
#include <conio.h>
#include <vector>
#include <complex>
#include <time.h>

using namespace std;

typedef complex<double> base;

void fft(vector<base>& a, bool invert) {
    int n = (int)a.size();

    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j >= bit; bit >>= 1)
            j -= bit;
        j += bit;
        if (i < j)
            swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * M_PI / len * (invert ? -1 : 1);
        base wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            base w(1);
            for (int j = 0; j < len / 2; ++j) {
                base u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if (invert)
        for (int i = 0; i < n; ++i)
            a[i] /= n;
}

long float* FIR_Filter(long float x[], long float y[], int size_n) {
    const int N = 20;
    long float fc = (1 + 0.4) / (20);  //50(50.0 + 20.0)/(4000.0)
    long float h[N] = { 0 };
    long float a = 0.5; // Параметр окна Гаусса

    long float hp = 2 * M_PI * fc;
    long float w = exp(-(2.0f * pow(N, 2.0f) / 4.0f) / (a * (float)N * a * (float)N));
    long float norm = 0;
    h[0] = hp * w;
    norm = h[0];

    for (int i = 1; i < N; i++) {
        hp = sinl(2 * M_PI * fc * i) / (M_PI * i);
        w = exp((2.0f * (i - pow(N, 2.0f) / 4.0f)) / (a * (float)N * a * (float)N)); // весовая функция для Гауссова окна
        h[i] = hp * w;
        norm += h[i];
        cout << hp << ' ' << w << endl;
    }
    /*
    norm = 0;
    for (int i = 0; i < N; i++) {
        if (i == 0) 
            hp = 2 * M_PI * fc;
        else 
            hp = sin(2 * M_PI * fc * i) / (M_PI * i);
        
        w = 0.42f - 0.5f * cosl((2.0 * M_PI * i) / (N - 1)) + 0.08f * cosl((4.0 * M_PI * i) / (N - 1)); // весовая функция Блекмена
        h[i] = hp * w;
        norm += hp;
    }
    */
    for (int i = 0; i < N; i++)
        h[i] /= norm;

    for (int i = 0; i < size_n; i++) {
        y[i] = 0.0;
        for (int j = 0; j < N - 1; j++) {
            if (i - j >= 0) {
                y[i] += h[j] * x[i - j];
            }
        }
        //out << y[i] << '\n';
    }
    return y;
}
