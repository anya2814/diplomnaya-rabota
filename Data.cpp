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

void Data::GetMoleculScatterCoef(std::vector<std::map<int, double>>& mol_koef_scat, double* waves)
{
    std::map<int, double> temp_map;
    std::ifstream R;
    std::string a;
    double H, Hprev, T, P, sr, sp_ip1[5], sp_i[5], value;
    R.open("Molekulyarnie parametri.txt");
    if (R.is_open()) {
        for (int i = 0; i < 3; i++) {
            R >> a;
        }

        // получаем первое значение для усреднения
        R >> H; R >> T; R >> P;
        for (int i = 0; i < 5; i++) {
            sr = 1 / (waves[i] * waves[i] * (938.076 * waves[i] * waves[i] - 10.8426));
            sp_ip1[i] = sr * P / 1013 * 273.16 / T;
        }

        Hprev = H;
        // получаем второе значение для усреднения, отдельно чтобы вставить векторы
        R >> H; R >> T; R >> P;
        for (int i = 0; i < 5; i++) {
            sr = 1 / (waves[i] * waves[i] * (938.076 * waves[i] * waves[i] - 10.8426));
            sp_i[i] = sp_ip1[i];
            sp_ip1[i] = sr * P / 1013 * 273.16 / T;
            value = 0.5 * (sp_i[i] + sp_ip1[i]);
            temp_map.insert(std::pair<int, double>(Hprev, value));
            mol_koef_scat.push_back(temp_map);
            temp_map.clear();
        }

        Hprev = H;
        for (int i = 3; i < 24; i++) {
            R >> H; R >> T; R >> P;
            for (int j = 0; j < 5; j++) {
                sr = 1 / (waves[j] * waves[j] * (938.076 * waves[j] * waves[j] - 10.8426));
                sp_i[j] = sp_ip1[j];
                sp_ip1[j] = sr * P / 1013 * 273.16 / T;
                value = 0.5 * (sp_i[j] + sp_ip1[j]);
                mol_koef_scat[j].insert(std::pair<int, double>(Hprev, value));
            }
            Hprev = H;
        }

        // последний слой где не делаем усреднение
        for (int j = 0; j < 5; j++) {
            sr = 1 / (waves[j] * waves[j] * (938.076 * waves[j] * waves[j] - 10.8426));
            sp_ip1[j] = sr * P / 1013 * 273.16 / T;
            mol_koef_scat[j].insert(std::pair<int, double>(Hprev, sp_ip1[j]));
        }
        R.close();
    }
}

void Data::getKoefOsl(std::vector<std::map<int, double>> &koef_osl, std::vector<std::map<int, double>> &alb_rass, double* waves)
{
    std::map<int, double> temp_map;
    int imass[5];
    int pos = 0;
    double read, h_next;
    std::ifstream H;
    H.open("AERO_MOD.txt");
    if (!H.is_open()) return;
    for (int i = 0; i < 27; i++) {
        H >> read;
        if (read == waves[pos]) {
            imass[pos] = i;
            pos++;
        }
    }
    H >> h_next;
    pos = 0;
    for (int i = 0; i < 27; i++) {
        H >> read;
        if (i == imass[pos]) {
            temp_map.insert(std::pair<int, double>(h_next, read));
            koef_osl.push_back(temp_map);
            temp_map.clear();
            pos++;
        }
    }
    for (int j = 1; j < 5; j++) {
        H >> h_next;
        pos = 0;
        for (int i = 0; i < 27; i++) {
            H >> read;
            if (i == imass[pos])
            {
                koef_osl[pos].insert(std::pair<int, double>(h_next, read));
                pos++;
            }
        }
    }

    H >> h_next;
    pos = 0;
    for (int i = 0; i < 27; i++) {
        H >> read;
        if (i == imass[pos]) {
            temp_map.insert(std::pair<int, double>(h_next, read));
            alb_rass.push_back(temp_map);
            temp_map.clear();
            pos++;
        }
    }
    for (int j = 1; j < 3; j++) {
        H >> h_next;
        pos = 0;
        for (int i = 0; i < 27; i++) {
            H >> read;
            if (i == imass[pos]) 
            {
                alb_rass[pos].insert(std::pair<int, double>(h_next, read));
                pos++;
            }
        }
    }
    H.close();
}

// вспомогательная функция для P6
// получение массива углов
float* Data::getM(float* angles)
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
            angles[i] = cos(read * PI / 180.);
            for (int j = 0; j < 5; j++) {
                H >> read;
            }
        }
        H.close();
    }

    return angles;
}
