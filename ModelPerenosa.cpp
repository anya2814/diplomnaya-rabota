#include "ModelPerenosa.h"

// получение случайного вещественного числа от 0 до 1
double ModelPerenosa::GetA() {
    double a;
    a = (rand() % 32767) / 32767.;
    return a;
}

double ModelPerenosa::sumUp = 0;
double ModelPerenosa::sumLow = 0;
int ModelPerenosa::L[5] = { 0, 0, 0, 0, 0 };

void ModelPerenosa::Lcount(int waveNum)
{
    L[waveNum]++;
}

double ModelPerenosa::GetSumUp()
{
    return sumUp;
}

int ModelPerenosa::GetL(int waveNum)
{
    return L[waveNum];
}

void ModelPerenosa::SetSum0()
{
    sumLow = 0;
    sumUp = 0;
}

double* ModelPerenosa::GetFi(double* fi) {
    double a1, a2;
    double w1 = 1, w2 = 1, d0 = 0;
    while (d0 == 0) {
        while ((w1 * w1 + w2 * w2) > 1)
        {
            a1 = GetA(); a2 = GetA();
            w1 = 1 - 2 * a1; w2 = 1 - 2 * a2;
        }
        d0 = w1 * w1 + w2 * w2;
        fi[0] = w1 / sqrt(d0);       // косинус фи
        fi[1] = w2 / sqrt(d0);       // синус фи
        w1 = 1; w2 = 1;
    }
        
    return fi;
}

// выбор начального распределения соответственно плотности распределения источника
/*double* ModelPerenosa::P1st_point(double* abc) {
    double* fi = new double[2];
    fi = GetFi(fi);
    abc[2] = 1 - 2 * GetA();
    double m = 1 - abc[2] * abc[2];

    abc[0] = fi[0];
    abc[1] = fi[1];

    return abc;
}*/

// выбор начального распределения соответственно плотности распределения источника
double* ModelPerenosa::P1st_point(double* abc) {
    double* fi = new double[2];
    fi = GetFi(fi);

    double cost = GetA();   // косинус тета
    double cosf = fi[0], sinf = fi[1]; // косинус и синус фи 
    double sint = sqrt(1 - cost * cost);    // синус тета
    
    abc[0] = sint*cosf;
    abc[1] = sint*sinf;
    abc[2] = cost;

    return abc;
}

// выбор длины свободного пробега l
void ModelPerenosa::P2length(int Lnum, double** d, double* xyz, double* abc, double* lopt, bool f) {
    double a = GetA();
    double l = -log(a);    //  
    double c = abc[2];          // косинус угла к поверхности Земли

    while ((a == 0) || (a == 1))
        a = GetA();
    l = -log(a);
    // длина полета больше чем оптическая длина поделить на косинус зенитного угла
    if ((l > (lopt[Lnum] / abc[2])) && (abc[2] > 0))
        Lcount(Lnum);

}

// функция для моделирования процесса переноса
void ModelPerenosa::ModPer(int Lnum, double** d, double* lopt) {
    double* abc = new double[3], * xyz = new double[3];

    for (int i = 0; i < 3; i++)
        abc[i] = 0;
    for (int i = 0; i < 3; i++)
        xyz[i] = 0;
    bool f = true;

    abc = P1st_point(abc);

    P2length(Lnum, d, xyz, abc, lopt, f);
}

int* ModelPerenosa::NModPer(int* t, int Lnum, double** d, double* lopt)
{
    for (int i = 0; i < kol; i++) {
        ModPer(Lnum, d, lopt);
    }

    return t;
}


void ModelPerenosa::Modelirovanie(double* waves, double** d, double* lopt)
{
    int* t = new int[3];
    for (int i = 0; i < 5; i++)
    {
        SetSum0();
        t = NModPer(t, i, d, lopt);
    }
    delete[]t;
}
