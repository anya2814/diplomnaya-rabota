#pragma once

#include <iostream>
#include <cmath>
#include <fstream>

// double abc[3] - ������ a, b � �
// double xyz[3] - ������ x, y, z
// double fi[2] - ������ cos ��, sin ��

static const int N = 204; // ����� �������� �������� F � m � �����
const double h = 30; // ������� ������� z
const int kol = 1000000; // ���������� ������������ �������� ������

class ModelPerenosa
{
    static double sumUp;
    static double sumLow;
    double GetA();                       // ��������� ���������� ����� � ��������� �� 0 �� 1 (����������� �������������)
    double getMa(float* mass, double** F, int Lnum, double a = ((rand() % 1001) / 1000.));
    double* GetFi(double* fi, double m = 1);            // ��������������� ������� ��� P1 � P7, ���������� �������� � ������ ��� ������ ��������� ����� � ��������� ��������� ����������� �������
    void CrossUp(double add);            // ���� ����������� ������� �������� � ����� 1/|(ns, w)|
    void CrossLow(double add);         // ���� ����������� ������ �������� � ����� 1/|(ns, w)|
    double* GetIzotr(double* abc);      // ����� ����������� ��� ����������� �������������
    double* GetLambert(double* abc);      // ����� ����������� ��� �������������� �������������
    double P2length(int Lnum, double** d, double z, double* abc, double pp);                   // ����� ����� ���������� ������� l
    double* P3P4calcul(double*, double*, double, double*);   // �������� ������ �� �����, ���������� ��������� ��������� ����� ������������
    bool P5type(int Lnum, double** d, double* xyz);                      // ����� ���� ������������ (���������� ��� ���������)
    double* P7napravl(float* mass, double** F, int Lnum, double* abc);      // �������� ��������� ����������� �������
    void Cout_xyz(double* xyz);
    void OutToFile(double** tBig, double* waves);

public:
    double GetSumUp();
    double GetSumLow();
    void SetSum0();
    void CountK(int* t);
    int ModPer(float* mass, double** F, int Lnum, double** d, double pp);
    int* NModPer(int* t, float* mass, double** F, int Lnum, double** d, double pp);
    void Modelirovanie(float* mass, double** F, double* waves, double** d, double pp);
};

