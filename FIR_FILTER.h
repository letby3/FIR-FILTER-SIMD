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
#include <string>

using namespace std;

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

class FirFilter {

private:
    float rm_maxx_end_fl[9][8] = { {0, 0, 0, 0, 0, 0, 0, 0},
                               {1, 0, 0, 0, 0, 0, 0, 0},
                               {1, 1, 0, 0, 0, 0, 0, 0},
                               {1, 1, 1, 0, 0, 0, 0, 0},
                               {1, 1, 1, 1, 0, 0, 0, 0},
                               {1, 1, 1, 1, 1, 0, 0, 0},
                               {1, 1, 1, 1, 1, 1, 0, 0},
                               {1, 1, 1, 1, 1, 1, 1, 0},
                               {1, 1, 1, 1, 1, 1, 1, 1} };

    __m256 rm_mass_end[9] = { _mm256_load_ps(&rm_maxx_end_fl[0][0]),
                              _mm256_load_ps(&rm_maxx_end_fl[1][0]),
                              _mm256_load_ps(&rm_maxx_end_fl[2][0]),
                              _mm256_load_ps(&rm_maxx_end_fl[3][0]),
                              _mm256_load_ps(&rm_maxx_end_fl[4][0]),
                              _mm256_load_ps(&rm_maxx_end_fl[5][0]),
                              _mm256_load_ps(&rm_maxx_end_fl[6][0]),
                              _mm256_load_ps(&rm_maxx_end_fl[7][0]),
                              _mm256_load_ps(&rm_maxx_end_fl[8][0]) };

public:

    typedef complex<double> base;

    struct dataset {
        float N = 128; //Точек отчеста входного сигнала 
        float fs = 10; //Частота дискретизации (в сек)
        float X_MAX = N / fs; //Длина сигнала в сек
        float* x; //Входной сигнал
        float* x_reverse;
        float* y; //Выходной сигнал
        int it = 0, it1 = N-1;

        __m256 N_simd;
        __m256 fs_simd;
        __m256 x_simd[126];
        __m256 y_simd[126];
    };
    
    ofstream stat_out;
    ofstream out;
    dataset main_config;

    FirFilter(string fileName, string stat_fileName) {
        main_config.x = new float[main_config.N];
        main_config.y = new float[main_config.N];
        main_config.x_reverse = new float[main_config.N];
        out.open(fileName);
        stat_out.open(stat_fileName);
    }

    FirFilter(float N, float fs, float X_MAX, float* x, float* y) {
        main_config.N = N;
        main_config.fs = fs;
        main_config.X_MAX = X_MAX;
        main_config.x = x;
        main_config.y = y;
    }

    void MakeNoise(float amp_min, float amp_max) {
        float* X_MAX = &main_config.X_MAX;
        float* fs = &main_config.fs;
        int it = main_config.it, it1 = main_config.it1;
        srand(time(NULL));
        for (float i = 0; i <= *X_MAX; i += 1.0 / *fs, it++, it1--) {
            //float X = ((float)(abs(rand() % 12)) / amp_min - amp_max) * cos(2 * M_PI * 0.5 * i); // Rand-ый sin-ый входной сигнал
            float X = (((float)rand() / (float)RAND_MAX) * (amp_max - amp_min) + amp_min) * cos(2 * M_PI * 0.5 * i);
            main_config.x[it] = X;
            main_config.x_reverse[it1] = X;         
        }
    }
    
    void OutXinGraph(string py_path, string GraphForFIR_path) {
        string python_x_param = py_path + " " + GraphForFIR_path + " graph1";
        char *python_x_param_char;
        for (int i = 0; i < main_config.N; i++) {
            python_x_param = python_x_param + " " + to_string(main_config.x[i]);
        }
        python_x_param_char = new char[python_x_param.size()];
        for (int i = 0; i < python_x_param.size(); i++) {
            python_x_param_char[i] = python_x_param[i];
        }
        system(python_x_param_char);
        delete python_x_param_char;
    }

    void OutYinGraph(string py_path, string GraphForFIR_path) {
        string python_y_param = py_path + " " + GraphForFIR_path + " graph2";
        char* python_y_param_char;
        for (int i = 0; i < main_config.N; i++) {
            python_y_param = python_y_param + " " + to_string(main_config.y[i]);            
        }
        python_y_param_char = new char[python_y_param.size()];
        for (int i = 0; i < python_y_param.size(); i++) {
            python_y_param_char[i] = python_y_param[i];
        }
        system(python_y_param_char);
        delete python_y_param_char;
    }

    void OutXYinGraph(string py_path, string GraphForFIR_path) {
        string python_xy_param = py_path + " " + GraphForFIR_path + " graph3";
        char* python_xy_param_char;
        for (int i = 0; i < main_config.N; i++) {
            python_xy_param = python_xy_param + " " + to_string(main_config.x[i]);
        }
        for (int i = 0; i < main_config.N; i++) {
            python_xy_param = python_xy_param + " " + to_string(main_config.y[i]);
        }
        python_xy_param_char = new char[python_xy_param.size()];
        
        for (int i = 0; i < python_xy_param.size(); i++) {
            python_xy_param_char[i] = python_xy_param[i];
        }
        system(python_xy_param_char);
        delete python_xy_param_char;
    }
    
    void start_second_test(string fileName, string stat_fileName) {
        out.open(fileName);
        stat_out.open(stat_fileName);
        srand(time(NULL));
        float N = main_config.N; //Точек отчеста входного сигнала 
        float fs = main_config.fs; //Частота дискретизации (в сек)
        float X_MAX = main_config.X_MAX; //Длина сигнала в сек
        float* x = main_config.x; //Входной сигнал
        float* x_reverse = main_config.x;
        float* y = main_config.y; //Выходной сигнал
        int it = main_config.it, it1 = main_config.it1;

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
            x_reverse[it1] = X;
            out << X << endl;           
        }
        for (int lenght_filter = 8; lenght_filter <= 2048; lenght_filter *= 2) {
            unsigned long long Fir_Filter_evg = 0,
                               Fir_Filter_SIMD_evg = 0;
        
            //FIR_Filter(x, y, N, lenght_filter);
        
            //Вывод сигналов от разных функций для сверки
            for (int i = 0; i < N; i++) {
                out << y[i] << endl; // Вывод выходного сигнала 
                y[i] = 0;
            }

            FIR_Filter_SIMD(x_reverse, y, N, lenght_filter);

            for (int i = 0; i < N; i++) {
                out << y[i] << endl; // Вывод выходного сигнала 
            }
        
            // проводим 1000 тестов для одного сигнала и замеряем время работы функции
            for (int num_test = 0; num_test < 10000;num_test++) {
                muTimer Fir_Filter;
                //FIR_Filter(x, y, N, lenght_filter);
                Fir_Filter_evg += Fir_Filter.stop().duration();
            
                muTimer Fir_Filter_SIMD;
                FIR_Filter_SIMD(x_reverse, y, N, lenght_filter);
                Fir_Filter_SIMD_evg +=  Fir_Filter_SIMD.stop().duration();
            
            }
            stat_out << (float)Fir_Filter_evg / 10000.0 << ' ' << (float)Fir_Filter_SIMD_evg / 10000.0 << ' '
                     << lenght_filter <<endl;
            cout << (float)Fir_Filter_evg / 10000.0 << ' ' << (float)Fir_Filter_SIMD_evg / 10000.0 << ' '
                << lenght_filter << endl;
        }
    }

    void start_first_test(string fileName) {       
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

        //FIR_Filter(x, y, N, 16); // КИХ фильтр нижних частот

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
        const int N_fl = inp_n; //Длина фильтра
        const int N_fl8 = N_fl / 8; //Длина фильтра /8
        float* hp = new float[N_fl]; // Импульсная характеристика
        float* w = new float[N_fl]; // Вес фильтра
        float fc_fl = (0.4 + 0.6) / 20.0; // аналогично 
        float a = 0.5; // Параметр окна Гаусса
    
    
        __m256* h = new __m256[N_fl8]; // Импульсная характеристика обр сигнала
        __m256* x_simd = new __m256[N_fl8]; //x_simd - часть входного сигнала в обратном порядке

    
        for (int i = 0; i < N_fl; i++) { // Расчёт веса и импульсной окна
            if (i == 0) {
                hp[i] = 2.0 * M_PI * fc_fl;
                w[i] = exp(-(2.0f * pow((float)N_fl, 2.0f) / 4.0f) / (a * (float)N_fl * a * (float)N_fl));
            }
            else {
                hp[i] = sinl(2 * M_PI * fc_fl * i) / (M_PI * i); //ФНЧ
                w[i] = exp((2.0f * (i - pow((float)N_fl, 2.0f) / 4.0f)) / (a * (float)N_fl * a * (float)N_fl)); // весовая функция для Гауссова окна
            }
        }
    
        __m256 Test = _mm256_setzero_ps();
        for (int i = 0, j = 0; i < N_fl; i += 8, j++) {
            h[j] = _mm256_mul_ps( _mm256_load_ps(&hp[i]), _mm256_load_ps(&w[i])); // Импульсная выходного. Расчёт
            Test = _mm256_add_ps(Test, h[j]); // Сумма импульсного 
        }
        float sum = 0;
        float* norm = (float*) &Test;

        for (int i = 0; i < 8; i++) {
            sum += norm[i]; // Сумма для нормирования
        }
        __m256 div = _mm256_set1_ps(sum);
        for (int i = 0; i < N_fl/8; i++) {
            h[i] = _mm256_div_ps(h[i], div); // Нормирование
        }

        __m256 sm = _mm256_setzero_ps(); // Сумма h[j] * x[i - j] для y[i]
        for (int i = 0; i < size_n; i ++) {

            sm = _mm256_setzero_ps();
            y[i] = 0;
            // -> вырезаем x[i] ... x[i + N_fl] (Для изначального сигнала) для реверсного x_r[size_n - i - 1] ... x_r[min(size_n - i - 1 + N_fl, size_n)] 
            for (int j = 0; j < N_fl/8; j++) {
                x_simd[j] = _mm256_setzero_ps();
                x_simd[j] = _mm256_load_ps(&x[size_n - i - 1 + j * 8]);
            }
            x_simd[min(i/8, N_fl/8 - 1)] = _mm256_mul_ps(x_simd[min(i / 8, N_fl / 8 - 1)], rm_mass_end[i%8 + 1]);
            //Домнажаем, чтобы избавиться от лишних значений в конце массива (Зануляем)
            // <-
        
            // -> FIR формула
            for (int j = 0; j <= min(i / 8, N_fl / 8 - 1); j++) {
                sm = _mm256_add_ps(sm, _mm256_mul_ps(h[j], x_simd[j]));
            }
            float* Sm = ((float*)&sm);
            float y_sm = Sm[0] + Sm[1] + Sm[2] + Sm[3] +
                         Sm[4] + Sm[5] + Sm[6] + Sm[7];
            y[i] = y_sm;
            // <-
        }

    
    }

    void FIR_Filter(int inp_n) {
        float *x = main_config.x;
        float* y = main_config.y; 
        int size_n = main_config.N;        
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
            h[i] /= norm; // Нормируем
        }

        for (int i = 0; i < size_n; i++) { //FIR filter
            y[i] = 0.0;
            for (int j = 0; j < N - 1; j++) {
                if (i - j >= 0) {
                    y[i] += h[j] * x[i - j];
                }
            }
        }
    }
};