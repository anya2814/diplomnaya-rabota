#pragma once

#include <iostream>
#include <cmath>

// double abc[3] - ������ a, b � �
// double xyz[3] - ������ x, y, z
// double fi[2] - ������ cos ��, sin ��

static const int N = 204; // ����� �������� �������� F � m � �����
const double h = 30; // ������� ������� z
const int kol = 100000; // ���������� ������������ �������� ������
#define PI 3.14159265

class ModelPerenosa
{
    static double sumUp;
    static double sumLow;
    static int L[5];
    double GetA();                       // ��������� ���������� ����� � ��������� �� 0 �� 1 (����������� �������������)
    double getMa(float* mass, double** F, int Lnum, double a = ((rand() % 1001) / 1000.));
    double* GetFi(double* fi, double m = 1);            // ��������������� ������� ��� P1 � P7, ���������� �������� � ������ ��� ������ ��������� ����� � ��������� ��������� ����������� �������
    void CrossUp(double add);            // ���� ����������� ������� �������� � ����� 1/|(ns, w)|
    void CrossLow(double add);         // ���� ����������� ������ �������� � ����� 1/|(ns, w)|
    void Lcount(int waveNum);          // ������� ������� ��� ����� ������ ���������� ������� ��������� ��� ������ �����
    double* P1st_point(double* abc);      // ����� ��������� ����� �������������� ��������� ������������� ���������
    double P2length(int Lnum, double** d, double* xyz, double* abc, double* lopt, bool f);                   // ����� ����� ���������� ������� l
    double* P3P4calcul(double*, double*, double, double);   // �������� ������ �� �����, ���������� ��������� ��������� ����� ������������
    bool P5type(int Lnum, double** d, double* xyz);                      // ����� ���� ������������ (���������� ��� ���������)
    double* P7napravl(double* abc, double m);      // �������� ��������� ����������� �������
    void Cout_xyz(double* xyz);

public:
    double GetSumUp();
    double GetSumLow();
    void SetSum0();
    int GetL(int waveNum);
    void CountK(int* t);
    int ModPer(float* mass, double** F, int Lnum, double** d, double* lopt);
    int* NModPer(int* t, float* mass, double** F, int Lnum, double** d, double* lopt);
    void Modelirovanie(float* mass, double** F, double* waves, double** d, double* lopt);
};

