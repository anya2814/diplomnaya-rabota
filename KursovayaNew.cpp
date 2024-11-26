// Kursovaya2903.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <time.h>
#include <map>
#include "ModelPerenosa.h"
#include "Data.h"

int main()
{
    srand(time(NULL));
    setlocale(LC_ALL, "rus");

    Data objData;
    double* waves = new double[5];
    std::vector<std::map<int, double>> koef_osl; // коэффициент ослабления потока
    std::vector<std::map<int, double>> alb_rass; // альбедо однократного рассеяния
    
    float* mass = new float[N];
    double** F = new double* [N];
    for (int i = 0; i < N; i++)
        F[i] = new double[5];
    waves = objData.getWaves(waves);    // длины волн
    mass = objData.getM(mass);      // массив углов
    F = objData.getF(F, mass);      // функция распределения угла рассеяния
    objData.getKoefOsl(koef_osl, alb_rass, waves);

    ModelPerenosa objModel;

    double pp = 0.35;   // задаем альбедо подстилающей поверхности

    // моделирование процессов переноса
    objModel.Modelirovanie(mass, F, waves, koef_osl, alb_rass, pp);

    delete[]mass;
    delete[]waves;
    for (int i = 0; i < N; i++)
        delete[]F[i];
    delete[]F;

    return 0;
}
