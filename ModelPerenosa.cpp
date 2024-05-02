#include "ModelPerenosa.h"

// получение случайного вещественного числа от 0 до 1
double ModelPerenosa::GetA() {
    double a;
    a = (rand() % 1001) / 1000.;
    return a;
}

double ModelPerenosa::sumUp = 0;
double ModelPerenosa::sumLow = 0;

// вспомогательная функция для P6
// выбор значения косинуса угла рассеяния m
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

// вспомогательная функция для P1 и P7
// нахождение косинуса и синуса для выбора начальной точки и пересчета координат направления пробега
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
    d0 = (w1 * w1 + w2 * w2) / m;
    fi[0] = w1 / sqrt(d0);     
    fi[1] = w2 / sqrt(d0);      
    w1 = 1; w2 = 1;

    return fi;
}

// учет пересечений верхней площадки с весом 1/|(ns, w)|
void ModelPerenosa::CrossUp(double add)
{
    sumUp = sumUp + 1.0 / (kol * abs(add));
}

double ModelPerenosa::GetSumUp()
{
    return sumUp;
}

// учет пересечений нижней площадки с весом 1/|(ns, w)|
void ModelPerenosa::CrossLow(double add)
{
    sumLow = sumLow + 1.0 / (kol * abs(add));
}

double ModelPerenosa::GetSumLow()
{
    return sumLow;
}

void ModelPerenosa::SetSum0()
{
    sumLow = 0;
    sumUp = 0;
}

// выбор начальной точки соответственно плотности распределения источника
double* ModelPerenosa::P1st_point(double* abc) {
    double* fi = new double[2];
    abc[2] = 0;
    while (
        (abc[2]==1) || (abc[2] == -1)) {
        abc[2] = 1 - 2 * GetA();
    }
    double m = 1 - abc[2] * abc[2];
    fi = GetFi(fi, m);

    abc[0] = fi[0];
    abc[1] = fi[1];

    return abc;
}

// выбор длины свободного пробега l
double ModelPerenosa::P2length(int Lnum, double** d, double z, double* abc) {
    int curr_ht = z / 1;        // слой в котором находится частица
    double a = GetA(), b = GetA(), l;
    double ln_prev = -log(a), ln_new, l_sum = 0;    //  
    double c = abc[2];          // косинус угла к поверхности Земли

    if (c == 0) return (ln_prev / d[Lnum][curr_ht]);  // если частица летит горизонтально

    if (c < 0) {
        // первый слой
        l = (z - curr_ht * 1.0) / abs(c);
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
            if (b <= 0.1) { c = -c; abc[2] = c; } 
            else return(-2); // произошло поглощение частицы поверхностью Земли
        }
    }

    else {
        // первый слой
        l = ((curr_ht * 1.0 + 1) - z) / c;
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

// проверка вылета из среды 
// вычисление координат очередной точки столкновения
double* ModelPerenosa::P3P4calcul(double* xyz, double* abc, double l, double c) {
    if (c != abc[2]) {
        xyz[0] = xyz[0] + abc[0] * l;
        xyz[1] = xyz[1] + abc[1] * l;
        l = l - xyz[2] / c;
        xyz[2] = abc[2] * l;
    }
    else {
        xyz[0] = xyz[0] + abc[0] * l;
        xyz[1] = xyz[1] + abc[1] * l;
        xyz[2] = xyz[2] + abc[2] * l;
    };

    return xyz;
}

// выбор типа столкновения (поглощение или рассеяние)
bool ModelPerenosa::P5type(int Lnum, double** d, double* xyz) {
    double a = GetA();
    int curr_ht = xyz[2] / 1;
    double Ps = d[Lnum][100 + curr_ht];// / d[Lnum][curr_ht];

    if (a < Ps) return 0; // произошло рассеяние
    else return 1; // произошло поглощение
}

// пересчет координат направления пробега
double* ModelPerenosa::P7napravl(double* abc, double m) {
    double* fi = new double[2], abct[3]; 
    double c = 0;
    for (int i = 0; i < 3; i++)
        abct[i] = abc[i];
    while ((c == 1) || (c == -1)) {
        fi = GetFi(fi);

        abc[0] = abct[0] * m - (abct[1] * fi[1] + abct[0] * abct[2] * fi[0]) * sqrt((1 - m * m) / (1 - abct[2] * abct[2]));
        abc[1] = abct[1] * m + (abct[0] * fi[1] - abct[1] * abct[2] * fi[0]) * sqrt((1 - m * m) / (1 - abct[2] * abct[2]));
        abc[2] = abct[2] * m + (1 - abct[2] * abct[2]) * fi[0] * sqrt((1 - m * m) / (1 - abct[2] * abct[2]));
        c = abc[2];
    }

    delete[]fi;
    return abc;
}

void ModelPerenosa::Cout_xyz(double* xyz) {

    std::cout << "x = " << xyz[0]
        << ", y = " << xyz[1]
        << ", z = " << xyz[2] << ";" << std::endl;
}

// функция для моделирования процесса переноса
int ModelPerenosa::ModPer(float* mass, double** F, int Lnum, double** d) {
    double* abc = new double[3], * xyz = new double[3];
    double l, m, c;

    for (int i = 0; i < 3; i++)
        abc[i] = 0;
    for (int i = 0; i < 3; i++)
        xyz[i] = 0;

    abc = P1st_point(abc);

    for (;;) {
        c = abc[2];
        l = P2length(Lnum, d, xyz[2], abc);

        if (l == -1)
        {
            // Произошел вылет за пределы среды через верхнюю границу
            if (abc[2] == 0)
                int p = 0;
            CrossUp(abc[2]);
            delete[]abc;
            delete[]xyz;
            return 1;
        }
        //cout << "l = " << l << endl;

        if (l == -2)
        {
            // Произошло поглощение частицы поверхностью Земли
            CrossLow(abc[2]);
            delete[]abc;
            delete[]xyz;
            return 0;
        }

        xyz = P3P4calcul(xyz, abc, l, c);
        //Cout_xyz(xyz);
        /*if (xyz[2] < 0)
        {
            xyz[2] = -xyz[2];
            abc[2] = abc[2] + 0.5;
            // Произошло отражение частицы от поверхности земли
        }*/
        /*if (xyz[2] > h)
        {
            // Произошел вылет за пределы среды через верхнюю границу
            CrossUp(abc[2]);
            delete[]abc;
            delete[]xyz;
            return 1;
        }*/
        if (P5type(Lnum, d, xyz)) {
            // Произошло поглощение
            delete[]abc;
            delete[]xyz;
            return 2;
        }
        m = getMa(mass, F, Lnum);
        abc = P7napravl(abc, m);

    }

}

int* ModelPerenosa::NModPer(int* t, float* mass, double** F, int Lnum, double** d)
{
    int k;
    for (int i = 0; i < 3; i++)
        t[i] = 0;
    for (int i = 0; i < kol; i++) {
        k = ModPer(mass, F, Lnum, d);
        t[k]++;
    }

    return t;
}

void ModelPerenosa::CountK(int* t)
{
    std::cout << "Коэффициент вылета через нижнюю границу: k1 = " << t[0] / (kol * 1.0) << std::endl;
    std::cout << "Коэффициент вылета через верхнюю границу: k2 = " << t[1] / (kol * 1.0) << std::endl;
    std::cout << "Коэффициент поглощений: k3 = " << t[2] / (kol * 1.0) << std::endl;
}

// вывод результатов в файл
void ModelPerenosa::OutToFile(double** tBig, double* waves)
{
    std::ofstream out;        
    out.open("Results.txt");      // открываем файл для записи
    if (out.is_open())
    {
        out << "Waves\tKtAbsorption\tKtUpCross\tKtUpCrossWeight\tKtLowCross\tKtLowCrossWeight" << std::endl;
        for (int i = 0; i < 5; i++) {
            out << waves[i] << '\t' << tBig[i][0] / (kol * 1.0) << '\t' << tBig[i][1] / (kol * 1.0) << '\t' << tBig[i][2] / (kol * 1.0);
            out << '\t' << tBig[i][3] / (kol * 1.0) << '\t' << tBig[i][4] / (kol * 1.0) << std::endl;
        }
    }
    out.close();
}

void ModelPerenosa::Modelirovanie(float* mass, double** F, double* waves, double** d)
{
    int* t = new int[3];
    double** tBig = new double* [5];
    for (int i = 0; i < 5; i++)
        tBig[i] = new double[5];
    for (int i = 0; i < 5; i++)
    {
        SetSum0();
        std::cout << "Данные для длины волны l=" << waves[i] << " мкм: " << std::endl << std::endl;
        t = NModPer(t, mass, F, i, d);
        std::cout << "Произошло " << t[0] << " поглощений частиц поверхностью Земли. " << std::endl;
        std::cout << "Произошло " << t[1] << " вылетов за пределы среды через верхнюю границу. " << std::endl;
        std::cout << "Произошло " << t[2] << " поглощений. " << std::endl;

        std::cout << std::endl;
        CountK(t);

        std::cout << std::endl;
        std::cout << "Коэффициент пересечения верхней границы с учетом веса: " << GetSumUp() << std::endl;
        std::cout << "Коэффициент поглощения нижней границей с учетом веса: " << GetSumLow() << std::endl << std::endl << std::endl;
        
        tBig[i][0] = t[2];
        tBig[i][1] = t[1];
        tBig[i][2] = GetSumUp();
        tBig[i][3] = t[0];
        tBig[i][4] = GetSumLow();
    }

    OutToFile(tBig, waves);

    delete[]t;
    for (int i = 0; i < 5; i++)
        delete[]tBig[i];
    delete[]tBig;
}
