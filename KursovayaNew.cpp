// Kursovaya2903.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <time.h>
#include "ModelPerenosa.h"
#include "Data.h"

int main()
{
    srand(time(NULL));
    setlocale(LC_ALL, "rus");

    Data objData;
    double* waves = new double[5];
    double** d = new double* [5];  // коэффициент ослабления потока - первые 100
    for (int i = 0; i < 5; i++) // альбедо однократного рассеяния - 100-200-ые
        d[i] = new double[200]; // делаю для каждой из 100 высот
    float* mass = new float[N];
    double** F = new double* [N];
    for (int i = 0; i < N; i++)
        F[i] = new double[5];
    waves = objData.getWaves(waves);
    mass = objData.getM(mass);
    F = objData.getF(F, mass);
    d = objData.getd(d, waves);

    ModelPerenosa objModel;

    // моделирование процессов переноса
    objModel.Modelirovanie(mass, F, waves, d);

    delete[]mass;
    delete[]waves;
    for (int i = 0; i < 2; i++)
        delete[]d[i];
    delete[]d;
    for (int i = 0; i < N; i++)
        delete[]F[i];
    delete[]F;

    return 0;
}
