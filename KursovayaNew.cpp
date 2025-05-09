// Kursovaya2903.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <windows.h>
#include <locale.h>
#include <time.h>
#include <map>
#include "ModelPerenosa.h"
#include "Data.h"

int main()
{
    srand(time(NULL));
    setlocale(LC_ALL, "Russian");
    SetConsoleOutputCP(1251);
    SetConsoleCP(1251);

    // переменные класса для вызова функций
    Data objData;
    ModelPerenosa objModel;

    // создаем массивы
    double* waves = new double[5];                // длины волн
    std::vector<std::map<int, double>> koef_osl; // коэффициент ослабления потока
    std::vector<std::map<int, double>> alb_rass; // альбедо однократного рассеяния
    std::vector<std::map<int, double>> mol_koef_rass; // альбедо однократного рассеяния
    
    // N - число заданных значений F и m в файле
    float* angles = new float[N];

    // заполняем массивы
    waves = objData.getWaves(waves);    // длины волн
    angles = objData.getM(angles);      // массив углов
    objData.GetMoleculScatterCoef(mol_koef_rass, waves);
    objData.getKoefOsl(koef_osl, alb_rass, waves);

    // задаем альбедо подстилающей поверхности
    double pp = 0;
    std::cout << "Vvedite znacheniye albedo podstilayushey poverhnosti:" << std::endl;
    std::cin >> pp;
    // задаем альбедо подстилающей поверхности
    // s - specular (зеркальное)
    // i - isotropic (изотропное)
    // l - lambertian (ламбертовское)
    int type = 0;
    std::cout << "Vyberite tip otrazheniya:" << std::endl << "1 - zerkalnoe" << std::endl << "2 - izotropnoe" << std::endl << "3 - lambertovskoe" << std::endl;
    std::cin >> type;

    // моделирование процессов переноса
    for (pp = 1; pp <= 1;) {
        objModel.Modelirovanie(angles, waves, mol_koef_rass, koef_osl, alb_rass, pp, type);
        pp = pp + 0.05;
    }

    // освобождение памяти
    delete[]angles;
    delete[]waves;

    return 0;
}
