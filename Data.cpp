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

void Data::getMoleculScatterCoef(std::vector<std::map<int, double>>& mol_koef_scat, double* waves)
{
    std::map<int, double> temp_map;
    std::ifstream R;
    std::string a;
    int H, Hprev;
    double T, P, sr, sp_ip1[5], sp_i[5], value;
    R.open("Molekulyarnie parametri.txt");
    if (R.is_open()) {
        for (int i = 0; i < 3; i++) {
            R >> a;
        }

        // получаем первое значение для усреднения
        R >> H; R >> T; R >> P;
        for (int i = 0; i < 5; i++) {
            sr = 1.0 / (waves[i] * waves[i] * (938.076 * waves[i] * waves[i] - 10.8426));
            sp_ip1[i] = sr * P / 1013 * 273.16 / T;
        }

        Hprev = H;
        // получаем второе значение для усреднения, отдельно чтобы вставить векторы
        R >> H; R >> T; R >> P;
        for (int i = 0; i < 5; i++) {
            sr = 1.0 / (waves[i] * waves[i] * (938.076 * waves[i] * waves[i] - 10.8426));
            sp_i[i] = sp_ip1[i];
            sp_ip1[i] = sr * P / 1013 * 273.16 / T;
            value = 0.5 * (sp_i[i] + sp_ip1[i]);
            temp_map.insert(std::pair<int, double>(Hprev, value));
            mol_koef_scat.push_back(temp_map);
            temp_map.clear();
        }

        Hprev = H;
        for (int i = 3; i < 25; i++) {
            R >> H; R >> T; R >> P;
            for (int j = 0; j < 5; j++) {
                sr = 1.0 / (waves[j] * waves[j] * (938.076 * waves[j] * waves[j] - 10.8426));
                sp_i[j] = sp_ip1[j];
                sp_ip1[j] = sr * P / 1013 * 273.16 / T;
                value = 0.5 * (sp_i[j] + sp_ip1[j]);
                mol_koef_scat[j][Hprev] = value;
            }
            Hprev = H;
        }

        // последний слой где не делаем усреднение
        for (int j = 0; j < 5; j++) {
            sr = 1.0 / (waves[j] * waves[j] * (938.076 * waves[j] * waves[j] - 10.8426));
            sp_ip1[j] = sr * P / 1013 * 273.16 / T;
            mol_koef_scat[j][Hprev] = sp_ip1[j];
        }
        R.close();
    }
}

void Data::getMoleculConsCoef(std::vector<std::map<int, double>>& mol_cons, double* waves)
{
    std::map<int, double> temp_map;
    int pos = 0;
    double read, h_next;

    std::ifstream file("mol_cons.txt");
    if (!file.is_open()) return;

    std::string header_line;
    std::getline(file, header_line);

    mol_cons.clear();
    mol_cons.resize(5);

    int H;
    double val0, val1, val2, val3, val4;

    // Читаем строки, пока файл не кончится
    while (file >> H >> val0 >> val1 >> val2 >> val3 >> val4) {
        mol_cons[0][H] = val0;
        mol_cons[1][H] = val1;
        mol_cons[2][H] = val2;
        mol_cons[3][H] = val3;
        mol_cons[4][H] = val4;
    }

    file.close();
}

#include <cmath>
#include <vector>
#include <map>
#include <fstream>

void Data::getKoefOsl(std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double* waves)
{
    std::ifstream in("AERO_MOD.txt");
    if (!in.is_open()) return;

    // 27 длин волн (в мкм) из файла
    const int M = 27;
    std::vector<double> w(M);
    for (int i = 0; i < M; ++i) in >> w[i];

    // Подготовим выходные структуры под 5 волн
    koef_osl.clear();
    alb_rass.clear();
    koef_osl.resize(5);
    alb_rass.resize(5);

    struct Bracket {
        int iL = 0;
        int iR = 0;
        double t = 0.0; // 0..1
    };

    Bracket br[5];

    // Для каждой из твоих 5 волн находим индексы слева/справа для интерполяции
    for (int j = 0; j < 5; ++j) {
        double lam = waves[j];

        int iL = 0, iR = 0;

        if (lam <= w[0]) {
            iL = iR = 0;
        }
        else if (lam >= w[M - 1]) {
            iL = iR = M - 1;
        }
        else {
            // ищем отрезок [w[i], w[i+1]], куда попадает lam
            for (int i = 0; i < M - 1; ++i) {
                if (w[i] <= lam && lam <= w[i + 1]) {
                    iL = i;
                    iR = i + 1;
                    break;
                }
            }
        }

        double t = (iL == iR) ? 0.0 : (lam - w[iL]) / (w[iR] - w[iL]);
        br[j] = { iL, iR, t };
    }

    // ---------- Блок EXTINCTION: 6 строк ----------
    // Формат: H + 27 чисел
    for (int row = 0; row < 6; ++row) {
        double H;
        if (!(in >> H)) return;

        std::vector<double> vals(M);
        for (int i = 0; i < M; ++i) in >> vals[i];

        int h_key = (int)std::lround(H);

        for (int j = 0; j < 5; ++j) {
            int L = br[j].iL;
            int R = br[j].iR;
            double t = br[j].t;

            double value = vals[L] + (vals[R] - vals[L]) * t;
            koef_osl[j][h_key] = value;
        }
    }

    // ---------- Блок SSA (альбедо однократного рассеяния): 4 строки ----------
    for (int row = 0; row < 4; ++row) {
        double H;
        if (!(in >> H)) return;

        std::vector<double> vals(M);
        for (int i = 0; i < M; ++i) in >> vals[i];

        int h_key = (int)std::lround(H);

        for (int j = 0; j < 5; ++j) {
            int L = br[j].iL;
            int R = br[j].iR;
            double t = br[j].t;

            double value = vals[L] + (vals[R] - vals[L]) * t;
            alb_rass[j][h_key] = value;
        }
    }
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
        for (int i = 0; i < 204; i++)
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

bool Data::readMolAbsKdistFile(const std::string& filename, KDistData& out)
{
    std::ifstream in(filename);
    if (!in.is_open()) return false;

    int Nexp = 0;
    double lam_left_nm = 0, lam_mid_nm = 0, lam_right_nm = 0;

    if (!(in >> Nexp >> lam_left_nm >> lam_mid_nm >> lam_right_nm)) return false;
    if (Nexp != 4) return false;

    out.lambda_um = lam_mid_nm / 1000.0; // nm -> um
    out.k_layers.clear();
    out.k_layers.resize(4);

    for (int i = 0; i < 4; ++i) {
        if (!(in >> out.Ci[i])) return false;
    }

    double z1 = 0, z2 = 0;
    double k[4]{ 0,0,0,0 };

    while (in >> z1 >> z2 >> k[0] >> k[1] >> k[2] >> k[3]) {
        int z_key = (int)std::lround(z1); // ключ = нижняя граница слоя
        for (int m = 0; m < 4; ++m) {
            out.k_layers[m][z_key] = k[m];
        }
    }

    return true;
}