#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "Data.h"

// double abc[3] - ������ a, b � �
// double xyz[3] - ������ x, y, z
// double fi[2] - ������ cos ��, sin ��

#define PI 3.14159265359

static const int N = 204; // ����� �������� �������� F � m � �����
const double h = 30; // ������� ������� z
const int kol = 500000; // ���������� ������������ �������� ������

struct IndicatrixTable {
    std::vector<double> lambda_um;                 // 5 значений из первой строки
    std::vector<double> theta_deg;                 // N=204
    std::vector<double> mu;                        // cos(theta)
    std::vector<std::vector<double>> P;            // P[j][i] : [5][N]
};

class ModelPerenosa
{
    static double sumUp;
    static double sumLow;
    
    std::map<double, KDistData> molAbsByLambda;  // все таблицы, ключ = lambda_um
    double tauMolPath[4]{ 0,0,0,0 };               // τ1..τ4 для одного фотона

    double GetA();                       // ��������� ���������� ����� � ��������� �� 0 �� 1 (����������� �������������)
    double getMa(double height, float* angles, double** F, double a = ((rand() % 1001) / 1000.));
    double* GetFi(double* fi, double m = 1);            // ��������������� ������� ��� P1 � P7, ���������� �������� � ������ ��� ������ ��������� ����� � ��������� ��������� ����������� �������
    void CrossUp(double* xyz, double* abc, double molConsWeight);            // ���� ����������� ������� �������� � ����� 1/|(ns, w)|
    void CrossLow(double* xyz, double* abc, double molConsWeight);         // ���� ����������� ������ �������� � ����� 1/|(ns, w)|
    void GetIzotr(double* abc);      // ����� ����������� ��� ����������� �������������
    void GetLambert(double* abc);      // ����� ����������� ��� �������������� �������������
    void findA(double *a, double e1_old[], double e2_old[], double e3_old[], double e_new[]);
    int P2length(int Lnum, double lambda_um, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& koef_osl, double* xyz, double* abc, double pp, int type);                   // ����� ����� ���������� ������� l
    bool Reflection(int Lnum, double* xyz, double* abc, double pp, int type = 1);
    std::pair<int, double> GetTequat(double* xyz, double* abc, double R); // ������� t �� ����������� ���������
    bool EarthReflType(double pp);
    bool P5type(int Lnum, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double* xyz);                      // ����� ���� ������������ (���������� ��� ���������)
    double* P7napravl(float* mass, double** F, int Lnum, double* abc, double height);      // �������� ��������� ����������� �������
    double** getF(double** F, float* angles, int Lnum, std::map<int, double>& mol_koef_rass, std::map<int, double>& koef_osl, std::map<int, double>& alb_rass);
    void Cout_xyz(double* xyz);
    void OutToFile(double** tBig, double* waves, double pp);
    double GetWeight(double* xyz, double* abc); // ����� ��� ��� �����������
    void AccumulateMolAbsTau(const KDistData* kdAbs, double ht_km, double ds_km);
    bool ReadIndicatrix(const std::string& filename, IndicatrixTable& tab);
    static double InterpByLambda(double lambda_um, const std::vector<double>& lam, const double p_at_lam[5]);

public:
    double GetSumUp();
    double GetSumLow();
    void SetSum0();
    void CountK(int* t);

    void AddMolAbsTable(const KDistData& kd) {
        molAbsByLambda[kd.lambda_um] = kd;
    }

    const KDistData* GetMolAbsForWave(double lambda_um) const;
    double MolAbsWeight(double lambda_um) const;

    int ModPer(float* mass, double lambda_um, double** F, int Lnum, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& mol_cons, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp, int type);
    int* NModPer(int* t, double lambda_um, float* mass, double** F, int Lnum, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& mol_cons, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp, int type);
    void Modelirovanie(float* mass, double* waves, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& mol_cons, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp, int type);
};

