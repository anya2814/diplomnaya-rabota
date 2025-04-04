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

    // переменные класса для вызова функций
    Data objData;
    ModelPerenosa objModel;

    // создаем массивы
    double* waves = new double[5];                // длины волн
    std::vector<std::map<int, double>> koef_osl; // коэффициент ослабления потока
    std::vector<std::map<int, double>> alb_rass; // альбедо однократного рассеяния
    
    // N - число заданных значений F и m в файле
    float* angles = new float[N];

    double** F = new double* [N];       // массив значений эмпирической функции распределения направлений рассеяния
    for (int i = 0; i < N; i++)
        F[i] = new double[5];

    // заполняем массивы
    waves = objData.getWaves(waves);    // длины волн
    angles = objData.getM(angles);      // массив углов
    F = objData.getF(F, angles);      // функция распределения угла рассеяния
    objData.getKoefOsl(koef_osl, alb_rass, waves);

    // задаем альбедо подстилающей поверхности
    double pp = 0;
    std::cout << "Введите значение альбедо подстилающей поверхности:" << std::endl;
    std::cin >> pp;
    // задаем альбедо подстилающей поверхности
    // s - specular (зеркальное)
    // i - isotropic (изотропное)
    // l - lambertian (ламбертовское)
    int type = 0;
    std::cout << "Выберите тип отражения:" << std::endl << "1 - зеркальное" << std::endl << "2 - изотропное" << std::endl << "3 - ламбертовское" << std::endl;
    std::cin >> type;

    // моделирование процессов переноса
    objModel.Modelirovanie(angles, F, waves, koef_osl, alb_rass, pp, type);

    // освобождение памяти
    delete[]angles;
    delete[]waves;
    for (int i = 0; i < N; i++)
        delete[]F[i];
    delete[]F;

    return 0;
}
