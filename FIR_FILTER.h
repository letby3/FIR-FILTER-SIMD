#include <iostream>
#include <cmath>
#include <fstream>
#include <conio.h>
#include <vector>
#include <complex>
#include <time.h>
#include <numeric>
#include <immintrin.h>
#include <chrono>
#include <ctime>

using namespace std;

typedef complex<double> base;

void FIR_Filter(float x[], float y[], int size_n, int inp_n);

void FIR_Filter_SIMD(float x[], float y[], int size_n, int inp_n);

void fft(vector<base>& a, bool invert);

void start_first_test(string fileName);

void start_second_test(string fileName);

float rm_maxx_end_fl[9][8] = { {0, 0, 0, 0, 0, 0, 0, 0},
                               {1, 0, 0, 0, 0, 0, 0, 0},
                               {1, 1, 0, 0, 0, 0, 0, 0}, 
                               {1, 1, 1, 0, 0, 0, 0, 0}, 
                               {1, 1, 1, 1, 0, 0, 0, 0}, 
                               {1, 1, 1, 1, 1, 0, 0, 0}, 
                               {1, 1, 1, 1, 1, 1, 0, 0}, 
                               {1, 1, 1, 1, 1, 1, 1, 0},
                               {1, 1, 1, 1, 1, 1, 1, 1}};

__m256 rm_mass_end[9] = { _mm256_load_ps(&rm_maxx_end_fl[0][0]),
                          _mm256_load_ps(&rm_maxx_end_fl[1][0]),
                          _mm256_load_ps(&rm_maxx_end_fl[2][0]),
                          _mm256_load_ps(&rm_maxx_end_fl[3][0]),
                          _mm256_load_ps(&rm_maxx_end_fl[4][0]),
                          _mm256_load_ps(&rm_maxx_end_fl[5][0]),
                          _mm256_load_ps(&rm_maxx_end_fl[6][0]),
                          _mm256_load_ps(&rm_maxx_end_fl[7][0]),
                          _mm256_load_ps(&rm_maxx_end_fl[8][0])};

class muTimer
{
    using Clock = std::chrono::high_resolution_clock;
    bool active = false;
    Clock::duration   duration_;
    Clock::time_point start_ = Clock::now(), stop_ = Clock::now();

    muTimer(const muTimer&) = delete;
    muTimer& operator=(const muTimer&) = delete;
public:
    using ns = std::chrono::nanoseconds;
    using mks = std::chrono::microseconds;
    using ms = std::chrono::milliseconds;
    muTimer() { reset(); start(); }
    ~muTimer() = default;
    muTimer& reset()
    {
        duration_ = std::chrono::nanoseconds(0);
        active = false;
        return *this;
    }
    muTimer& start()
    {
        if (!active)
        {
            start_ = Clock::now();
            active = true;
        }
        return *this;
    }
    muTimer& stop()
    {
        if (active)
        {
            stop_ = Clock::now();
            duration_ += stop_ - start_;
            active = false;
        }
        return *this;
    }
    template<typename T = mks>
    unsigned long long duration()
    {
        return static_cast<unsigned long long>
            (std::chrono::duration_cast<T>(stop_ - start_).count());
    }
};

void start_second_test(string fileName) {
    
    ofstream out;
    out.open(fileName);
    srand(time(NULL));
    float N = 128; //Точек отчеста входного сигнала 
    float fs = 10; //Частота дискретизации (в сек)
    float X_MAX = N / fs; //Длина сигнала в сек
    float* x = new float[1001]; //Входной сигнал
    float* x_reverse = new float[1001];
    vector<float> xx;
    float* y = new float[1001]; //Выходной сигнал
    int it = 0, it1 = 127;

    __m256 N_simd;
    __m256 fs_simd;
    __m256 x_simd[126];
    __m256 y_simd[126];

    for (long float i = 0, j = 0; i <= X_MAX - 1.0/fs; i += 1.0 / fs, j++) {
        out << i << endl; //Вывод координаты ВХ сигнала в сек
    }

    for (float i = 0; i <= X_MAX; i += 1.0 / fs, it++, it1--) {
        float X = ((float)(abs(rand() % 12)) / 100.0 - 0.01) * cos(2 * M_PI * 0.5 * i); // Rand-ый sin-ый входной сигнал
        x[it] = X;
        xx.push_back(X);
        x_reverse[it1] = X;
        out << X << endl;
    }

    muTimer mt;
    FIR_Filter(x, y, N, 64);
    cout << mt.stop().duration() << " mks\n";

    for (int i = 0; i < N; i++) {
        out << y[i] << endl; // Вывод выходного сигнала 
        y[i] = 0;
    }

    muTimer mt1;
    FIR_Filter_SIMD(x_reverse, y, N, 64);
    cout << mt1.stop().duration() << " mks\n";

    for (int i = 0; i < N; i++) {
        out << y[i] << endl; // Вывод выходного сигнала 
    }
}

void start_first_test(string fileName) {
    ofstream out;
    out.open(fileName);
    srand(time(NULL));
    float N = 128; //Точек отчеста входного сигнала 
    float fs = 10; //Частота дискретизации (в сек)
    float X_MAX = N / fs - 1 / fs; //Длина сигнала в сек
    float* x = new float[1001]; //Входной сигнал
    float* y = new float[1001]; //Выходной сигнал
    int it = 0; //верный друг, верный счётчик
    vector<base> a, b; // ВХОД и ВЫХОД Фурье

    for (float i = 0, j = 0; i <= X_MAX; i += 1.0 / fs, j++) {
        out << i << endl; //Вывод координаты ВХ сигнала в сек
    }

    for (float i = 0; i <= X_MAX; i += 1.0 / fs, it++) {
        float X = ((long float)(abs(rand() % 12)) / 100.0 - 0.01) * cos(2 * M_PI * 0.5 * i); // Rand-ый sin-ый входной сигнал
        out << X << endl;
        x[it] = X;
        a.push_back(X);
    }

    fft(a, 0); //Фурье входного сигнала (Complex)

    FIR_Filter(x, y, N, 16); // КИХ фильтр нижних частот

    for (int i = 0; i < N; i++) {
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
}

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

void FIR_Filter_SIMD(float x[], float y[], int size_n, int inp_n) {
    const int N_fl = inp_n;
    const int N_fl8 = N_fl / 8;
    float* hp = new float[N_fl];
    float* w = new float[N_fl];
    float fc_fl = (0.4 + 0.6) / 20.0;
    float a = 0.5;
    __m256 N = _mm256_set1_ps((float)N_fl);
    __m256 fc = _mm256_set1_ps(fc_fl);
    __m256* h = new __m256[N_fl8];
    __m256* x_simd = new __m256[N_fl8];
    //__m256 a = _mm256_set1_ps(0.5);
    for (int i = 0; i < N_fl; i++) {
        if (i == 0) {
            hp[i] = 2.0 * M_PI * fc_fl;
            w[i] = exp(-(2.0f * pow((float)N_fl, 2.0f) / 4.0f) / (a * (float)N_fl * a * (float)N_fl));
        }
        else {
            hp[i] = sinl(2 * M_PI * fc_fl * i) / (M_PI * i); //ФНЧ
            w[i] = exp((2.0f * (i - pow((float)N_fl, 2.0f) / 4.0f)) / (a * (float)N_fl * a * (float)N_fl)); // весовая функция для Гауссова окна
        }
    }
    //__m256 Test = _mm256_mul_ps(_mm256_load_ps(&hp[8]), _mm256_load_ps(&w[8]));
    __m256 Test = _mm256_setzero_ps();
    for (int i = 0, j = 0; i < N_fl; i += 8, j++) {
        h[j] = _mm256_mul_ps( _mm256_load_ps(&hp[i]), _mm256_load_ps(&w[i]));
        Test = _mm256_add_ps(Test, h[j]);
    }
    float sum = 0;
    float* norm = (float*) &Test;
    for (int i = 0; i < 8; i++) {
        sum += norm[i];
    }
    __m256 div = _mm256_set1_ps(sum);
    for (int i = 0; i < N_fl/8; i++) {
        h[i] = _mm256_div_ps(h[i], div);
    }

    __m256 sm = _mm256_setzero_ps();
    for (int i = 0; i < size_n; i ++) {
        sm = _mm256_setzero_ps();
        y[i] = 0;
        for (int j = 0; j < N_fl/8; j++) {
            x_simd[j] = _mm256_setzero_ps();
            x_simd[j] = _mm256_load_ps(&x[size_n - i - 1 + j * 8]);
        }
        x_simd[min(i/8, N_fl/8 - 1)] = _mm256_mul_ps(x_simd[min(i / 8, N_fl / 8 - 1)], rm_mass_end[i%8 + 1]);
        
        for (int j = 0; j <= min(i / 8, N_fl / 8 - 1); j++) {
            sm = _mm256_add_ps(sm, _mm256_mul_ps(h[j], x_simd[j]));
        }
        float* Sm = ((float*)&sm);
        float y_sm = Sm[0] + Sm[1] + Sm[2] + Sm[3] +
                     Sm[4] + Sm[5] + Sm[6] + Sm[7];
        //cout << y_sm << endl;
        y[i] = y_sm;
    }
    
}

void FIR_Filter(float x[], float y[], int size_n, int inp_n) {
    const int N = inp_n; // Порядок фильтра
    float fc = (0.4 + 0.6) / (20);  // fc = (fp + fst)/(2*fs) := (полоса пропускания + п. заграждения)/(2 * частота дискретизации)
    float* h = new float[N]; //Импульсная характеристика
    float a = 0.5; // Параметр окна Гаусса

    float hp = 2 * M_PI * fc; // hp импульсная характеристика фильтра
    float w = exp(-(2.0f * pow(N, 2.0f) / 4.0f) / (a * (float)N * a * (float)N)); // w вес фильтра
    float norm = 0; // Для нормирования импульсной характеристики
    h[0] = hp * w;
    norm = h[0];

    for (int i = 1; i < N; i++) {
        hp = sinl(2.0 * M_PI * fc * (float)i) / (M_PI * (float)i); //ФНЧ
        w = exp((2.0f * ((float)i - pow(N, 2.0f) / 4.0f)) / (a * (float)N * a * (float)N)); // весовая функция для Гауссова окна
        h[i] = hp * w;
        norm += h[i];
        //cout << hp << ' ' << w << endl;
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
    for (int i = 0; i < N; i++) {
        h[i] /= norm;
    }

    for (int i = 0; i < size_n; i++) { //FIR filter
        y[i] = 0.0;
        for (int j = 0; j < N - 1; j++) {
            if (i - j >= 0) {
                //cout << i << ' ' << j << endl;
                y[i] += h[j] * x[i - j];
            }
        }
    }
}
