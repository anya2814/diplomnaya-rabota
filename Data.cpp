#include "Data.h"

// вспомогательная функция для P6
// чтение из файла длин волн
double* Data::getWaves(double* wavesLength)
{
    std::ifstream H;
    H.open("Indikatrisa.txt");
    if (H.is_open()) {
        for (int i = 0; i < 5; i++) {
            H >> wavesLength[i];
        }
        H.close();
    }

    return wavesLength;
}

double** Data::getd(double** d, double* waves)
{
    int* imass = new int[5];
    int pos = 0;
    double read, h_prev = 0, h_next;
    std::ifstream H;
    H.open("AERO_MOD.txt");
    if (!H.is_open()) return nullptr;
    for (int i = 0; i < 27; i++) {
        H >> read;
        if (read == waves[pos]) {
            imass[pos] = i;
            pos++;
        }
    }
    for (int j = 0; j < 30; j++) {
        H >> h_next;
        pos = 0;
        for (int i = 0; i < 27; i++) {
            H >> read;
            if (i == imass[pos]) {
                for (int k = h_prev; k < h_next; k++) {
                    d[pos][k] = read;
                }
                pos++;
            }
        }
        h_prev = h_next;
    }
    h_prev = 0;
    for (int j = 0; j < 30; j++) {
        H >> h_next;
        pos = 0;
        for (int i = 0; i < 27; i++) {
            H >> read;
            if (i == imass[pos]) {
                for (int k = h_prev; k < h_next; k++) {
                    d[pos][100+k] = read;
                }
                pos++;
            }
        }
        h_prev = h_next;
    }

    H.close();

    return d;
}

// вспомогательная функция для P6
// получение массива углов
float* Data::getM(float* mass)
{
    std::ifstream H;
    float read;
    H.open("Indikatrisa.txt");
    if (H.is_open()) {
        for (int i = 0; i < 5; i++) {
            H >> read;
        }
        for (int i = 0; i < N; i++)
        {
            H >> read;
            mass[i] = cos(read * PI / 180.);
            for (int j = 0; j < 5; j++) {
                H >> read;
            }
        }
        H.close();
    }

    return mass;
}

// вспомогательная функция для P6
// вычисление F по формуле трапеций
double** Data::getF(double** F, float* mass)
{
    double** ind = new double* [N];
    for (int i = 0; i < N; i++)
        ind[i] = new double[5];
    double sum[5];
    for (int i = 0; i < 5; i++)
        sum[i] = 0;
    std::ifstream H;
    float read;
    H.open("Indikatrisa.txt");
    if (H.is_open()) {
        for (int i = 0; i < 5; i++) {
            H >> read;
        }
        bool flag = false;
        for (int i = 0; i < N; i++)
        {
            H >> read;
            for (int j = 0; j < 5; j++) {
                H >> ind[i][j];
                if (flag)
                    sum[j] = sum[j] + (ind[i][j] + ind[i - 1][j]) / 2.0 * abs(mass[i] - mass[i - 1]);
                F[i][j] = sum[j];
            }
            flag = true;
        }
        H.close();
        for (int i = 0; i < N; i++)
            for (int j = 0; j < 5; j++)
                F[i][j] = F[i][j] / F[N - 1][j];
    }
    for (int i = 0; i < N; i++)
        delete[]ind[i];
    delete[]ind;

    //writeF(F);

    return F;
}

/*void Data::writeF(double** F)
{
    std::ofstream H;
    H.open("F.txt");
    H.clear();
    if (H.is_open())
    {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < 4; j++)
                H << F[i][j] << " ";
            H << F[i][4];
            H << std::endl;
        }
    }
    H.close();
}*/