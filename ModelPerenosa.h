#pragma once

#include <iostream>
#include <cmath>

// double abc[3] - массив a, b и с
// double xyz[3] - массив x, y, z
// double fi[2] - массив cos фи, sin фи

static const int N = 204; // число заданных значений F и m в файле
const double h = 30; // верхняя граница z
const int kol = 100000; // количество моделируемых пробегов частиц

class ModelPerenosa
{
    static double sumUp;
    static double sumLow;
    static int L[5];
    double GetA();                       // получение случайного числа в интервале от 0 до 1 (равномерное распределение)
    void Lcount(int waveNum);          // подсчет сколько раз длина больше оптической толщины атмосферы для данной длины
    double* P1st_point(double* abc);      // выбор начальной точки соответственно плотности распределения источника
    void P2length(int Lnum, double** d, double* xyz, double* abc, double* lopt, bool f);                   // выбор длины свободного пробега l

public:
    double GetSumUp();
    void SetSum0();
    int GetL(int waveNum);
    double* GetFi(double* fi);
    void ModPer(int Lnum, double** d, double* lopt);
    int* NModPer(int* t, int Lnum, double** d, double* lopt);
    void Modelirovanie(double* waves, double** d, double* lopt);
};

