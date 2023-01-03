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

    typedef complex<float> base;

    struct dataset {
        float N = 128; //Точек отчеста входного сигнала 
        float fs = 1000; //Частота дискретизации (в сек)
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
    ifstream in_read_wav;
    dataset main_config;

    FirFilter() {
        main_config.x = new float[main_config.N];
        main_config.y = new float[main_config.N];
        main_config.x_reverse = new float[main_config.N];
    }

    FirFilter(float N, float fs, float X_MAX, float* x, float* y) {
        main_config.N = N;
        main_config.fs = fs;
        main_config.X_MAX = X_MAX;
        main_config.x = x;
        main_config.y = y;
    }

    void UploadWavSignal(string py_path, string file_name) {
        string python_x_param = py_path + " " + file_name + " " + "readWav ";
        char* python_x_param_char;
        python_x_param_char = new char[python_x_param.size()];
        for (int i = 0; i < python_x_param.size(); i++) {
            python_x_param_char[i] = python_x_param[i];
        }
        system(python_x_param_char);
        delete python_x_param_char;
        in_read_wav.open("C:\\Users\\Пользователь\\source\\repos\\RpojectOOOSTC1\\RpojectOOOSTC1\\file_read.txt");
        in_read_wav >> main_config.fs;
        in_read_wav >> main_config.N;     
        main_config.X_MAX = main_config.N / main_config.fs;
        main_config.x = new float[main_config.N];
        main_config.y = new float[main_config.N];
        main_config.x_reverse = new float[main_config.N];
        for (int i = 0; i < main_config.N; i++)
            in_read_wav >> main_config.x[i];
    }

    void ExportWavSignalY(string py_path, string GraphForFIR_path, string file_name) {
        string python_x_param = py_path + " " + GraphForFIR_path + " writeWav " + file_name;
        char* python_x_param_char;
        out.close();
        out.open(file_name + ".txt");
        vector<int> wav_write_int16;
        float max_wav_y = 0;
        for (int i = 0; i < main_config.N; i++)
            max_wav_y = max(main_config.y[i], max_wav_y);
        for (int i = 0; i < main_config.N; i++)
            out << (int)main_config.y[i] << endl; //((int)((main_config.y[i] / max_wav_y) * 32767.0))

        python_x_param += " ";
        python_x_param_char = new char[python_x_param.size()];
        for (int i = 0; i < python_x_param.size(); i++) {
            python_x_param_char[i] = python_x_param[i];
        }
        system(python_x_param_char);
        delete python_x_param_char;
    }

    void MakeNoise(float amp_min, float amp_max, int sin_on_off) {
        float* X_MAX = &main_config.X_MAX;
        float* fs = &main_config.fs;
        int it = main_config.it, it1 = main_config.it1;
        srand(time(NULL));
        for (float i = 0; i <= *X_MAX; i += 1.0 / *fs, it++, it1--) {
            //float X = ((float)(abs(rand() % 12)) / amp_min - amp_max) * cos(2 * M_PI * 0.5 * i); // Rand-ый sin-ый входной сигнал
            float X = (((float)rand() / (float)RAND_MAX) * (amp_max - amp_min) + amp_min) * cos(2 * M_PI * 0.5 * i * sin_on_off);
            cout << X << endl;
            main_config.x[it] = X;
            main_config.x_reverse[it1] = X;         
        }
    }
    
    void OutXinGraph(string py_path, string GraphForFIR_path, string file_out) {
        string python_x_param = py_path + " " + GraphForFIR_path + " graph1 " + file_out + " ";
        char *python_x_param_char;
        out.close();
        out.open(file_out);        
        for (int i = 0; i < main_config.N; i++) {
            out << main_config.x[i] << endl;
        }        
        python_x_param_char = new char[python_x_param.size()];
        for (int i = 0; i < python_x_param.size(); i++) {
            python_x_param_char[i] = '\0';
            python_x_param_char[i] = python_x_param[i];                
        }                
        system(python_x_param_char);
        delete python_x_param_char;        
    }

    void OutYinGraph(string py_path, string GraphForFIR_path, string file_out) {
        string python_y_param = py_path + " " + GraphForFIR_path + " graph2 " + file_out + " ";
        char* python_y_param_char;
        out.close();
        out.open(file_out);
        for (int i = 0; i < main_config.N; i++) {
            out << (main_config.y[i]) << endl;
        }
        python_y_param_char = new char[python_y_param.size()];
        for (int i = 0; i < python_y_param.size(); i++) {
            python_y_param_char[i] = python_y_param[i];
        }
        system(python_y_param_char);
        delete python_y_param_char;
    }

    void OutXYinGraph(string py_path, string GraphForFIR_path, string file_out) {
        string python_xy_param = py_path + " " + GraphForFIR_path + " graph3 " + file_out + " ";
        char* python_xy_param_char;
        out.close();
        out.open(file_out);          
        for (int i = 0; i < main_config.N; i++) {
            out << main_config.x[i] << endl;
        }
        for (int i = 0; i < main_config.N; i++) {
            out << main_config.y[i] << endl;
        }
        python_xy_param_char = new char[python_xy_param.size()];
        
        for (int i = 0; i < python_xy_param.size(); i++) {
            python_xy_param_char[i] = python_xy_param[i];
        }
        system(python_xy_param_char);
        delete python_xy_param_char;
    }
    

    void FFT(vector<base>& a, bool invert) {       
        long long n = (int)a.size();

        for (long long i = 1, j = 0; i < n; ++i) {
            long long bit = n >> 1;
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

    void OutFrequencyResponseGraph(string py_path, string GraphForFIR_path, string file_out) {
        string python_frr_param = py_path + " " + GraphForFIR_path + " graph4 " + file_out + " ";
        char* python_frr_param_char;
        vector<base> x, y;
        out.close();
        out.open(file_out);
        int new_n = 1;
        while (1) {
            if (new_n * 2 > main_config.N)
                break;
            new_n *= 2;
        }        
        
        for (int i = 0; i < new_n;i++) {
            x.push_back(main_config.x[i]);
            y.push_back(main_config.y[i]);
        }
        FFT(x, 0);
        FFT(y, 0);        
        for (int i = 2; i < x.size(); i++) {            
            out << (sqrt(x[i].real() * x[i].real() + x[i].imag() * x[i].imag())) << endl;
        }
        for (int i = 2; i < y.size(); i++) {
            out << (sqrt(y[i].real() * y[i].real() + y[i].imag() * y[i].imag())) << endl;
        }        
        python_frr_param += to_string(x.size() + y.size()) + " ";
        python_frr_param_char = new char[python_frr_param.size()];
        for (int i = 0; i < python_frr_param.size(); i++) {
            python_frr_param_char[i] = python_frr_param[i];
        }
        system(python_frr_param_char);
        delete(python_frr_param_char);
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

    void FIR_Filter(int inp_n, float fp, float fst) {
        float *x = main_config.x;
        float* y = main_config.y; 
        int size_n = main_config.N;        
        const int N = inp_n; // Порядок фильтра
        float fc = (fp + fst) / (2 * main_config.fs);  // fc = (fp + fst)/(2*fs) := (полоса пропускания + п. заграждения)/(2 * частота дискретизации)
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