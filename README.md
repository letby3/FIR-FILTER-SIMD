# OOOSTC_ProblemsForApprentice
All Problems u can find in "задачи.docs"

2.2	Написать реализацию КИХ фильтра, оптимизированную с использованием 
расширений SIMD процессора x86_64. Сравнить производительность 
оптимизированной и неоптимизированной реализации в зависимости от длины 
импульсной характеристики фильтра.

Литература: 
1. https://habr.com/ru/post/460445/
2. https://hub.exponenta.ru/post/osnovy-tsifrovoy-obrabotki-signalov-achkh-i-fchkh-tsifrovye-filtry-kikh-i-bikh-filtry612
3. https://ppt-online.org/47360
4. http://www.dsplib.ru/
5. https://www.geokniga.org/bookfiles/geokniga-davydovlekciipocifrovojobrabotkesignalov03-effektgibbsaivesovie.pdf
6. https://wiki5.ru/wiki/Window_function
7. https://ru.stackoverflow.com/questions/789161/%D0%98%D0%B7%D0%BC%D0%B5%D1%80%D0%B5%D0%BD%D0%B8%D0%B5-%D0%B2%D1%80%D0%B5%D0%BC%D0%B5%D0%BD%D0%B8-%D0%B2%D1%8B%D0%BF%D0%BE%D0%BB%D0%BD%D0%B5%D0%BD%D0%B8%D1%8F-%D0%BF%D1%80%D0%BE%D0%B3%D1%80%D0%B0%D0%BC%D0%BC%D1%8B-%D0%B2-%D0%BD%D0%B0%D0%BD%D0%BE%D1%81%D0%B5%D0%BA%D1%83%D0%BD%D0%B4%D0%B0%D1%85
8. https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html

Решение:

Hardware: i5-10300H

OS: Ubuntu 18.04

Compiler: MK VS C++ 14, Python 3.10.8
(Графики для удобства строил в Python)

Файлы: "FiniteImpulseResponse.cpp", "FIR_FILTER.h", "GraphForFIR.py", "stat_out.txt", "graph1.txt"
![Figure_1](https://user-images.githubusercontent.com/22713938/204394177-108faf75-d434-439f-9eff-4dbd81b5fc6e.png)

  1.1. Графики с lenght_filter = 8, 16, 32, 64
На графиках с маленькой длиной фильтра наблюдается большой разброс 
(ссылаюсь на разные типы переменных (баг не пофиксил пока))
С большей длинной фильтра такого не наблюдается

![Figure_2](https://user-images.githubusercontent.com/22713938/204394439-5b42e661-57a7-416e-ad4e-2bf501c51cc7.png)

  1.2. Графики с lenght_filter = 128, 256, 512, 1024

![Figure_4](https://user-images.githubusercontent.com/22713938/204394456-f5be31e6-f135-4118-8317-59fd14b076c1.png)

  1.3 Графики синий-FIR_Filter, оранжевый-FIR_Filter_SIMD
