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

double* Data::lOpt(double* l, double** d, double* waves) {
    for (int i = 0; i < 5; i++)
        l[i] = 0;
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 100; j++)
            l[i] = l[i] + d[i][j];
    saveLOpt(waves, l);
    return l;
}

void Data::saveLOpt(double* waves, double* l) {
    std::ofstream out;
    out.open("Optical depth.txt");      // открываем файл для записи
    if (out.is_open())
    {
        out.clear();
        out << "Waves" << '\t' << "Optical depth" << std::endl;
        for (int i = 0; i < 5; i++)
            out << waves[i] << '\t' << l[i] << std::endl;
    }
    out.close();
}