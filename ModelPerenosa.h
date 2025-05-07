#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <fstream>

// double abc[3] - ������ a, b � �
// double xyz[3] - ������ x, y, z
// double fi[2] - ������ cos ��, sin ��

#define PI 3.14159265

static const int N = 204; // ����� �������� �������� F � m � �����
const double h = 30; // ������� ������� z
const int kol = 100000; // ���������� ������������ �������� ������

class ModelPerenosa
{
    static double sumUp;
    static double sumLow;
    double GetA();                       // ��������� ���������� ����� � ��������� �� 0 �� 1 (����������� �������������)
    double getMa(float* mass, double** F, int Lnum, double a = ((rand() % 1001) / 1000.));
    double* GetFi(double* fi, double m = 1);            // ��������������� ������� ��� P1 � P7, ���������� �������� � ������ ��� ������ ��������� ����� � ��������� ��������� ����������� �������
    void CrossUp(double* xyz, double* abc);            // ���� ����������� ������� �������� � ����� 1/|(ns, w)|
    void CrossLow(double* xyz, double* abc);         // ���� ����������� ������ �������� � ����� 1/|(ns, w)|
    void GetIzotr(double* abc);      // ����� ����������� ��� ����������� �������������
    void GetLambert(double* abc);      // ����� ����������� ��� �������������� �������������
    void findA(double *a, double e1_old[], double e2_old[], double e3_old[], double e_new[]);
    int P2length(int Lnum, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double* xyz, double* abc, double pp, int type);                   // ����� ����� ���������� ������� l
    bool Reflection(int Lnum, double* xyz, double* abc, double pp, int type = 1);
    std::pair<int, double> GetTequat(double* xyz, double* abc, double R); // ������� t �� ����������� ���������
    bool EarthReflType(double pp);
    bool P5type(int Lnum, std::vector<std::map<int, double>>& alb_rass, double* xyz);                      // ����� ���� ������������ (���������� ��� ���������)
    double* P7napravl(float* mass, double** F, int Lnum, double* abc);      // �������� ��������� ����������� �������
    double** getF(double** F, float* angles, int Lnum, std::map<int, double>& mol_koef_rass, std::map<int, double>& koef_osl, std::map<int, double>& alb_rass);
    void Cout_xyz(double* xyz);
    void OutToFile(double** tBig, double* waves, double pp);
    double GetWeight(double* xyz, double* abc); // ����� ��� ��� �����������

public:
    double GetSumUp();
    double GetSumLow();
    void SetSum0();
    void CountK(int* t);
    int ModPer(float* mass, double** F, int Lnum, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp, int type);
    int* NModPer(int* t, float* mass, double** F, int Lnum, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp, int type);
    void Modelirovanie(float* mass, double* waves, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp, int type);
};

