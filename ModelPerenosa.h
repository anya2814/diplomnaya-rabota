#pragma once

#include <iostream>
#include <cmath>

// double abc[3] - ������ a, b � �
// double xyz[3] - ������ x, y, z
// double fi[2] - ������ cos ��, sin ��

static const int N = 204; // ����� �������� �������� F � m � �����
const double h = 30; // ������� ������� z
const int kol = 100000; // ���������� ������������ �������� ������

class ModelPerenosa
{
    static double sumUp;
    static double sumLow;
    static int L[5];
    double GetA();                       // ��������� ���������� ����� � ��������� �� 0 �� 1 (����������� �������������)
    void Lcount(int waveNum);          // ������� ������� ��� ����� ������ ���������� ������� ��������� ��� ������ �����
    double* P1st_point(double* abc);      // ����� ��������� ����� �������������� ��������� ������������� ���������
    void P2length(int Lnum, double** d, double* xyz, double* abc, double* lopt, bool f);                   // ����� ����� ���������� ������� l

public:
    double GetSumUp();
    void SetSum0();
    int GetL(int waveNum);
    double* GetFi(double* fi);
    void ModPer(int Lnum, double** d, double* lopt);
    int* NModPer(int* t, int Lnum, double** d, double* lopt);
    void Modelirovanie(double* waves, double** d, double* lopt);
};

