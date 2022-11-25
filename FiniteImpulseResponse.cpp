#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <conio.h>

using namespace std;

ofstream out;
long float answ[2000];

void FIR_Filter(long float x[], long float y[], const int size_n);

int main(){
    out.open("graph.txt");
    long float N = 1000;
    long float X_MAX = 10;
    long float *x = new long float[2000];
    long float *y = new long float[2000];
    for (long float i = -X_MAX; i < X_MAX; i += (2 * X_MAX / N)) {
        out << i << endl;
    }
    for (long float i = -X_MAX; i < X_MAX; i += (2 * X_MAX / N)) {
        out << 100*sin(i) << endl;
    }

    FIR_Filter(x, y, 2*N);

    system("PAUSE");
    return 0;
}

void FIR_Filter(long float x[], long float y[], int size_n) {
    const int N = 20;
    long float fc = (50.0 + 20.0)/(4000.0);  //50
    long float h[N] = { 0 };
    long float a = 0.5; // Параметр окна Гаусса

    long float hp = 2 * M_PI * fc;
    long float w = exp(-(2.0f * pow(N, 2.0f) / 4.0f) / (a * (float)N * a * (float)N));
    long float norm = 0;
    h[0] = hp * w;
    norm = h[0];

    for (int i = 1; i < N; i++) {
        hp = sin(2 * M_PI * fc * i) / (M_PI * i);
        w = exp(-(2.0f * (i - pow(N, 2.0f) / 4.0f)) / (a * (float)N * a * (float)N));
        h[i] = hp * w;
        norm += h[i];
        //cout << hp << ' ' << w << endl;
    }

    for (int i = 0; i < N; i++)
        h[i] /= norm;

    for (int i = 0; i < size_n; i++) {
        y[i] = 0;
        for (int j = 0; j < N - 1; j++) {
            if (i - j >= 0) {
                y[i] += h[j] * x[i - j];
            }
        }
        out << y[i] << '\n';
    }
}
