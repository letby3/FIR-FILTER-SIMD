#define _USE_MATH_DEFINES
#include <iostream>
#include <conio.h>
#include "FIR_FILTER.h"
#include <immintrin.h>

using namespace std;


int main(){
    string py_path = "python";
    string FirFilter_path = "C:\\Users\\Пользователь\\source\\repos\\RpojectOOOSTC1\\OOOSTC_ProblemsForApprentice\\GraphForFIR.py";
    FirFilter test1 = FirFilter("graph1.txt", "stat_out.txt");
    test1.UploadWavSignal(py_path, FirFilter_path);
    test1.OutXinGraph(py_path, FirFilter_path, "C:\\Users\\Пользователь\\source\\repos\\RpojectOOOSTC1\\RpojectOOOSTC1\\graphX.txt");
    test1.FIR_Filter(1024, 9000, 7966);
    test1.OutFrequencyResponseGraph(py_path, FirFilter_path, "C:\\Users\\Пользователь\\source\\repos\\RpojectOOOSTC1\\RpojectOOOSTC1\\graphXYFR.txt");
    test1.ExportWavSignalY(py_path, FirFilter_path, "C:\\Users\\Пользователь\\source\\repos\\RpojectOOOSTC1\\RpojectOOOSTC1\\file_write");
    system("PAUSE");
    return 0;
}
