#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <fstream>

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
    double GetA();                       // получение случайного числа в интервале от 0 до 1 (равномерное распределение)
    double getMa(float* mass, double** F, int Lnum, double a = ((rand() % 1001) / 1000.));
    double* GetFi(double* fi, double m = 1);            // вспомогательная функция для P1 и P7, нахождение косинуса и синуса для выбора начальной точки и пересчета координат направления пробега
    void CrossUp(double add);            // учет пересечений верхней площадки с весом 1/|(ns, w)|
    void CrossLow(double add);         // учет пересечений нижней площадки с весом 1/|(ns, w)|
    void GetIzotr(double* abc, double *xyz);      // выбор направления для изотропного распределения
    void GetLambert(double* abc, double* xyz);      // выбор направления для ламбертовского распределения
    int P2length(int Lnum, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double* xyz, double* abc, double pp);                   // выбор длины свободного пробега l
    bool Reflection(double* xyz, double* abc);
    double GetTequat(double* xyz, double* abc, double R); // находим t из квадратного уравнения
    bool P5type(int Lnum, std::vector<std::map<int, double>>& alb_rass, double* xyz);                      // выбор типа столкновения (поглощение или рассеяние)
    double* P7napravl(float* mass, double** F, int Lnum, double* abc);      // пересчет координат направления пробега
    void Cout_xyz(double* xyz);
    void OutToFile(double** tBig, double* waves);

public:
    double GetSumUp();
    double GetSumLow();
    void SetSum0();
    void CountK(int* t);
    int ModPer(float* mass, double** F, int Lnum, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp);
    int* NModPer(int* t, float* mass, double** F, int Lnum, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp);
    void Modelirovanie(float* mass, double** F, double* waves, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp);
};

