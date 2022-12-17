#define _USE_MATH_DEFINES
#include <iostream>
#include <conio.h>
#include "FIR_FILTER.h"
#include <immintrin.h>

using namespace std;


int main(){
    string py_path = "C:\\Users\\letby\\OneDrive\\Work\\Projects\\GraphFIR\\venv\\Scripts\\python.exe";
    string FirFilter_path = "C:\\Users\\letby\\source\\repos\\FirFilterNewEdit\\FirFilterNewEdit\\GraphForFIR.py";
    FirFilter test1 = FirFilter("graph1.txt", "stat_out.txt");
    //test1.start_second_test("graph1.txt", "stat_out.txt");
    test1.MakeNoise(200.0, 20.0); 
    //test1.OutXinGraph(py_path, FirFilter_path);
    test1.FIR_Filter(32);
    //test1.OutYinGraph(py_path, FirFilter_path);
    test1.OutXYinGraph(py_path, FirFilter_path);
    system("PAUSE");
    return 0;
}
