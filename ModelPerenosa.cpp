#include "ModelPerenosa.h"

// ��������� ���������� ������������� ����� �� 0 �� 1
double ModelPerenosa::GetA() {
    double a;
    a = (rand() % 1001) / 1000.;
    return a;
}

double ModelPerenosa::sumUp = 0;
double ModelPerenosa::sumLow = 0;
int ModelPerenosa::L[5] = { 0, 0, 0, 0, 0 };

// ��������������� ������� ��� P6
// ����� �������� �������� ���� ��������� m
double ModelPerenosa::getMa(float* mass, double** F, int Lnum, double a)
{
    int leftI, rightI;
    double cosm;
    if (a == 0) cosm = mass[0];
    else for (int i = 1; i < N; i++)
        if (a <= F[i][Lnum]) {
            rightI = i - 1; leftI = i; i = N;
            cosm = mass[rightI] - (mass[rightI] * 1.0 - mass[leftI] * 1.0) * (a - F[rightI][Lnum]) / (F[leftI][Lnum] - F[rightI][Lnum]);
        }
    return cosm;
}

// ��������������� ������� ��� P1 � P7
// ���������� �������� � ������ ��� ������ ��������� ����� � ��������� ��������� ����������� �������
double* ModelPerenosa::GetFi(double* fi, double m) {
    double a1, a2;
    double w1 = 1, w2 = 1, d0;
    fi[1] = 1;
    if (m == 0) {
        fi[0] = 0;
        fi[1] = 0;
        return fi;
    }
    while ((w1 * w1 + w2 * w2) > m)
    {
        a1 = GetA(); a2 = GetA();
        w1 = 1 - 2 * a1; w2 = 1 - 2 * a2;
    }
    d0 = (w1 * w1 + w2 * w2)/m;
    fi[0] = w1 / sqrt(d0);       // �������
    fi[1] = w2 / sqrt(d0);       // �����
    w1 = 1; w2 = 1;
    
    if (abs((fi[0] * fi[0] + fi[1] * fi[1]) - m) > 0.0000001)
        int i=0;
    return fi;
}

void ModelPerenosa::Lcount(int waveNum)
{
    L[waveNum]++;
}

// ���� ����������� ������� �������� � ����� 1/|(ns, w)|
void ModelPerenosa::CrossUp(double add)
{
    sumUp = sumUp + 1.0 / (kol * abs(add));
}

double ModelPerenosa::GetSumUp()
{
    return sumUp;
}

// ���� ����������� ������ �������� � ����� 1/|(ns, w)|
void ModelPerenosa::CrossLow(double add)
{
    sumLow = sumLow + 1.0 / (kol * abs(add));
}

double ModelPerenosa::GetSumLow()
{
    return sumLow;
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

// ����� ��������� ����� �������������� ��������� ������������� ���������
double* ModelPerenosa::P1st_point(double* abc) {
    double a1, a2, a3;
    double w1 = 1, w2 = 1, w3 = 1, d0;

    while ((w1 * w1 + w2 * w2 + w3*w3) > 1)
    {
        a1 = GetA(); a2 = GetA(); a3 = GetA();
        w1 = 1 - 2*a1; w2 = 1 - 2 * a2; w3 = 1 - 2*a3;
    }
    d0 = w1 * w1 + w2 * w2 + w3 * w3;
    abc[0] = w1 / sqrt(d0);   
    abc[1] = w2 / sqrt(d0);       
    abc[2] = w3 / sqrt(d0);

    return abc;
}

/*double* ModelPerenosa::P1st_point(double* abc) {
    double* fi = new double[2];

    abc[2] = 1 - 2*GetA();
    double m = 1 - abc[2] * abc[2];
    fi = GetFi(fi, m);

    abc[1] = fi[0];
    abc[2] = fi[1];

    return abc;
}*/

// ����� ����� ���������� ������� l
double ModelPerenosa::P2length(int Lnum, double** d, double* xyz, double* abc, double* lopt, bool f) {
    int curr_ht = xyz[2] / 1;        // ���� � ������� ��������� �������
    double a = GetA(), l;
    double ln_prev = -log(a), ln_new, l_sum = 0;    //  
    double c = abc[2];          // ������� ���� � ����������� �����

    if (abs(abc[2] - 1) > 0.00001)
        int i = 0;
    if ((ln_prev > (lopt[Lnum] / c)) && f && (c > 0))
        Lcount(Lnum);

    if (c == 0) return (ln_prev / d[Lnum][curr_ht]);  // ���� ������� ����� �������������

    if (c < 0) {
        // ������ ����
        l = (xyz[2] - curr_ht * 1.0) / abs(c);
        ln_new = ln_prev - l * d[Lnum][curr_ht];

        if (ln_new <= 0) return (ln_prev / d[Lnum][curr_ht]);
        
        if (curr_ht == 0) { c = -c; abc[2] = c; }
        else {
            curr_ht--;
            for (; curr_ht >= 0;) {
                ln_prev = ln_new;
                l = 1 / abs(c);
                ln_new = ln_prev - l * d[Lnum][curr_ht];
                if (ln_new <= 0) return (l_sum + ln_prev / d[Lnum][curr_ht]);
                l_sum += l;
                curr_ht--;
            }
            c = -c; abc[2] = c;
        }
    }

    else {
        // ������ ����
        l = ((curr_ht * 1.0 + 1) - xyz[2]) / c;
        ln_new = ln_prev - l * d[Lnum][curr_ht];
        if (ln_new <= 0) return (ln_prev / d[Lnum][curr_ht]);
        l_sum += l;
        curr_ht++;
    }

    for (; curr_ht < 100;) {
        ln_prev = ln_new;
        l = 1 / c;
        if (d[Lnum][curr_ht] < 0.00000001)
            return -1;
        ln_new = ln_prev - l * d[Lnum][curr_ht];
        if (ln_new <= 0) return (l_sum + ln_prev / d[Lnum][curr_ht]);
        l_sum += l;
        curr_ht++;
    }
}

// �������� ������ �� ����� 
// ���������� ��������� ��������� ����� ������������
double* ModelPerenosa::P3P4calcul(double* xyz, double* abc, double l, double c) {
    if (c != abc[2]) {
        l = l - xyz[2] / c;
    }
    else {
        xyz[0] = xyz[0] + abc[0] * l;
        xyz[1] = xyz[1] + abc[1] * l;
        xyz[2] = xyz[2] + abc[2] * l;
    };

    return xyz;
}

// ����� ���� ������������ (���������� ��� ���������)
bool ModelPerenosa::P5type(int Lnum, double** d, double* xyz) {
    double a = GetA();
    int curr_ht = xyz[2] / 1;
    double Ps = d[Lnum][100 + curr_ht];// / d[Lnum][curr_ht];

    if (a < Ps) return 0; // ��������� ���������
    else return 1; // ��������� ����������
}

// �������� ��������� ����������� �������
double* ModelPerenosa::P7napravl(double* abc, double m) {
    double* fi = new double[2], abct[3]; 
    for (int i = 0; i < 3; i++)
        abct[i] = abc[i];
    fi = GetFi(fi);

    abc[0] = abct[0] * m - (abct[1] * fi[1] + abct[0] * abct[2] * fi[0]) * sqrt((1 - m * m) / (1 - abct[2] * abct[2]));
    abc[1] = abct[1] * m + (abct[0] * fi[1] - abct[1] * abct[2] * fi[0]) * sqrt((1 - m * m) / (1 - abct[2] * abct[2]));
    abc[2] = abct[2] * m + (1 - abct[2] * abct[2]) * fi[0] * sqrt((1 - m * m) / (1 - abct[2] * abct[2]));

    delete[]fi;
    return abc;
}

void ModelPerenosa::Cout_xyz(double* xyz) {

    std::cout << "x = " << xyz[0]
        << ", y = " << xyz[1]
        << ", z = " << xyz[2] << ";" << std::endl;
}

// ������� ��� ������������� �������� ��������
int ModelPerenosa::ModPer(float* mass, double** F, int Lnum, double** d, double* lopt) {
    double* abc = new double[3], * xyz = new double[3];
    double l, m, c;

    for (int i = 0; i < 3; i++)
        abc[i] = 0;
    for (int i = 0; i < 3; i++)
        xyz[i] = 0;
    bool f = true;

    abc = P1st_point(abc);

    for (;;) {

        c = abc[2];
        l = P2length(Lnum, d, xyz, abc, lopt, f);
        return 1;
        f = false;
        if (l == -1)
        {
            // ��������� ����� �� ������� ����� ����� ������� �������
            CrossUp(abc[2]);
            delete[]abc;
            delete[]xyz;
            return 1;
        }
        //cout << "l = " << l << endl;

        xyz = P3P4calcul(xyz, abc, l, c);
        //Cout_xyz(xyz);
        /*if (xyz[2] < 0)
        {
            xyz[2] = -xyz[2];
            abc[2] = abc[2] + 0.5;
            // ��������� ��������� ������� �� ����������� �����
        }*/
        /*if (xyz[2] > h)
        {
            // ��������� ����� �� ������� ����� ����� ������� �������
            CrossUp(abc[2]);
            delete[]abc;
            delete[]xyz;
            return 1;
        }*/
        if (P5type(Lnum, d, xyz)) {
            // ��������� ����������
            delete[]abc;
            delete[]xyz;
            return 2;
        }

        m = getMa(mass, F, Lnum);
        abc = P7napravl(abc, m);
    }

}

int* ModelPerenosa::NModPer(int* t, float* mass, double** F, int Lnum, double** d, double* lopt)
{
    int k;
    for (int i = 0; i < 3; i++)
        t[i] = 0;
    for (int i = 0; i < kol; i++) {
        k = ModPer(mass, F, Lnum, d, lopt);
        t[k]++;
    }

    return t;
}

void ModelPerenosa::CountK(int* t)
{
    std::cout << "����������� ������ ����� ������ �������: k1 = " << t[0] / (kol * 1.0) << std::endl;
    std::cout << "����������� ������ ����� ������� �������: k2 = " << t[1] / (kol * 1.0) << std::endl;
    std::cout << "����������� ����������: k3 = " << t[2] / (kol * 1.0) << std::endl;
}

void ModelPerenosa::Modelirovanie(float* mass, double** F, double* waves, double** d, double* lopt)
{
    int* t = new int[3];
    for (int i = 0; i < 5; i++)
    {
        SetSum0();
        std::cout << "������ ��� ����� ����� l=" << waves[i] << " ���: " << std::endl << std::endl;
        t = NModPer(t, mass, F, i, d, lopt);
        std::cout << "��������� " << t[0] << " ������� �� ������� ����� ����� ������ �������. " << std::endl;
        std::cout << "��������� " << t[1] << " ������� �� ������� ����� ����� ������� �������. " << std::endl;
        std::cout << "��������� " << t[2] << " ����������. " << std::endl;

        std::cout << std::endl;
        CountK(t);

        std::cout << std::endl;
        std::cout << "����������� ����������� ������� ������� � ������ ����: " << GetSumUp() << std::endl;
        std::cout << "����������� ����������� ������ ������� � ������ ����: " << GetSumLow() << std::endl << std::endl << std::endl;
    }
    delete[]t;
}
