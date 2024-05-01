#pragma once

#include <iostream>
#include <cmath>

// double abc[3] - массив a, b и с
// double xyz[3] - массив x, y, z
// double fi[2] - массив cos фи, sin фи

static const int N = 204; // число заданных значений F и m в файле
const double h = 30; // верхняя граница z
const int kol = 100000; // количество моделируемых пробегов частиц
#define PI 3.14159265

class ModelPerenosa
{
    static double sumUp;
    static double sumLow;
    static int L[5];
    double GetA();                       // получение случайного числа в интервале от 0 до 1 (равномерное распределение)
    double getMa(float* mass, double** F, int Lnum, double a = ((rand() % 1001) / 1000.));
    double* GetFi(double* fi, double m = 1);            // вспомогательная функция для P1 и P7, нахождение косинуса и синуса для выбора начальной точки и пересчета координат направления пробега
    void CrossUp(double add);            // учет пересечений верхней площадки с весом 1/|(ns, w)|
    void CrossLow(double add);         // учет пересечений нижней площадки с весом 1/|(ns, w)|
    void Lcount(int waveNum);          // подсчет сколько раз длина больше оптической толщины атмосферы для данной длины
    double* P1st_point(double* abc);      // выбор начальной точки соответственно плотности распределения источника
    double P2length(int Lnum, double** d, double* xyz, double* abc, double* lopt, bool f);                   // выбор длины свободного пробега l
    double* P3P4calcul(double*, double*, double, double);   // проверка вылета из среды, вычисление координат очередной точки столкновения
    bool P5type(int Lnum, double** d, double* xyz);                      // выбор типа столкновения (поглощение или рассеяние)
    double* P7napravl(double* abc, double m);      // пересчет координат направления пробега
    void Cout_xyz(double* xyz);

public:
    double GetSumUp();
    double GetSumLow();
    void SetSum0();
    int GetL(int waveNum);
    void CountK(int* t);
    int ModPer(float* mass, double** F, int Lnum, double** d, double* lopt);
    int* NModPer(int* t, float* mass, double** F, int Lnum, double** d, double* lopt);
    void Modelirovanie(float* mass, double** F, double* waves, double** d, double* lopt);
};

