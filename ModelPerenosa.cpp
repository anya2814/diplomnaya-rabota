#include "ModelPerenosa.h"

// получение случайного вещественного числа от 0 до 1
double ModelPerenosa::GetA() {
    double a;
    a = (rand() % 32767) / 32767.;
    return a;
}

double ModelPerenosa::sumUp = 0;
double ModelPerenosa::sumLow = 0;

// вспомогательная функция для P6
// выбор значения косинуса угла рассеяния m
double ModelPerenosa::getMa(float* angles, double** F, int Lnum, double a)
{
    int leftI, rightI;
    double cosm;
    if (a == 0) cosm = angles[0];
    else for (int i = 1; i < N; i++)
        if (a <= F[i][Lnum]) {
            rightI = i - 1; leftI = i; i = N;
            cosm = angles[rightI] - (angles[rightI] * 1.0 - angles[leftI] * 1.0) * (a - F[rightI][Lnum]) / (F[leftI][Lnum] - F[rightI][Lnum]);
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

// выбор направления для ламбертовского распределения
void ModelPerenosa::GetLambert(double* abc) {
    double* fi = new double[2];
    fi = GetFi(fi);

    double cost = 1;   // косинус тета
    while (cost == 1)
        cost = sqrt(GetA());

    double cosf = fi[0], sinf = fi[1]; // косинус фи
    double sint = sqrt(1 - cost * cost);    // синус тета

    abc[0] = sint * cosf;
    abc[1] = sint * sinf;
    abc[2] = cost;

    delete[]fi;
}

// выбор направления для изотропного распределения
void ModelPerenosa::GetIzotr(double* abc) {
    double* fi = new double[2];
    fi = GetFi(fi);

    double cost = 1;   // косинус тета
    while (cost == 1)
        cost = GetA();

    double cosf = fi[0], sinf = fi[1]; // косинус фи
    double sint = sqrt(1 - cost * cost);    // синус тета
    
    abc[0] = sint * cosf;
    abc[1] = sint * sinf;
    abc[2] = cost;

    delete[]fi;
}

// ПРЕДЫДУЩАЯ ФУНКЦИЯ выбор длины свободного пробега l
/*double ModelPerenosa::P2length(int Lnum, double** d, double z, double* abc, double pp) {
    int curr_ht = z / 1;        // слой в котором находится частица
    double a = GetA(), b = GetA(), l;
    double ln_prev = -log(a), ln_new, l_sum = 0;    //  
    double c = abc[2];          // косинус угла к поверхности Земли

    while ((a == 0) || (a == 1)) {
        a = GetA();
    }
    ln_prev = -log(a);

    if (c == 0) return (ln_prev / d[Lnum][curr_ht]);  // если частица летит горизонтально

    if (c < 0) {
        // первый слой
        l = (z - curr_ht * 1.0) / abs(c);
        ln_new = ln_prev - l * d[Lnum][curr_ht];
        l_sum += l;

        if (ln_new <= 0) return (ln_prev / d[Lnum][curr_ht]);

        if (curr_ht == 0) {
            if (b <= pp) {
                abc = GetLambert(abc);
                c = abc[2];
                if (c == 0) return (l_sum + ln_prev / d[Lnum][curr_ht]);  // если частица летит горизонтально
            }
            else return(-2); // произошло поглощение частицы поверхностью Земли
        }
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
            if (b <= pp) { 
                abc = GetLambert(abc);
                c = abc[2]; 
                if (c == 0) return (l_sum + ln_prev / d[Lnum][curr_ht]);  // если частица летит горизонтально
            }
            else return(-2); // произошло поглощение частицы поверхностью Земли
        }
    }

    else {
        // первый слой
        l = ((curr_ht * 1.0 + 1) - z) / c;
        if (d[Lnum][curr_ht] < 0.00000001)
            return -1;
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
}*/



// НОВАЯ ФУНКЦИЯ выбор длины свободного пробега l + проверка вылета из среды 
// вычисление координат очередной точки столкновения
int ModelPerenosa::P2length(int Lnum, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double* xyz, double* abc, double pp, int type) {
    double ht = sqrt(xyz[2] * xyz[2] + xyz[1] * xyz[1] + xyz[0] * xyz[0]) - 6371;
    double j = sqrt(xyz[2] * xyz[2] + xyz[1] * xyz[1] + xyz[0] * xyz[0]) - 6371;

    std::map<int, double>::iterator curr_ht_ko; // слой в котором находится частица для коэффициента ослабления  
                                                // (0 - 0-3 км, 1 - 3-13 км, 2 - 13-25 км, 3 - 25-35 км, 4 - 35-100 км, где ко = 0)
    for (std::map<int, double>::iterator it = koef_osl[Lnum].begin(); it != koef_osl[Lnum].end(); it++) {
        if (ht >= it->first) curr_ht_ko = it; // указатель на std::pair где находится текущий коэффициент ослабления    
    }

    double a = 0, R, l_opt;//, l_real = 0;    
    double c = abc[2];          // косинус угла к поверхности Земли
    std::pair<int, double> t = std::make_pair(-1, 1);            // количество корней и нужный корень

    while ((a == 0) || (a == 1)) {  // чтобы не брать логарифмы от 0 и 1
        a = GetA();
    }
    l_opt = -log(a);

    double temp = 0; // временная переменная для проверки

    t.second = 0;
    // цикл пока частица летит внутрь
    while (t.first == 2 || t.first == -1) {
        R = curr_ht_ko->first + 6371.0; // радиус сферы
        t = GetTequat(xyz, abc, R);
        if (alb_rass[Lnum].size() != 3)
            int k = 0;
        if (t.first == 2 && t.second != -1) {
            temp = l_opt - t.second * curr_ht_ko->second; // от общей оптической длины отнимаем сколько частица пролетает до слоя
            if (temp <= 0)
                t.second = l_opt / curr_ht_ko->second;
            else l_opt = l_opt - t.second * curr_ht_ko->second;
            //l_real = l_real + t;    // реальная длина пробега
            if ((sqrt(xyz[2] * xyz[2] + xyz[1] * xyz[1] + xyz[0] * xyz[0]) - 6371) < -0.01)
                j = sqrt(xyz[2] * xyz[2] + xyz[1] * xyz[1] + xyz[0] * xyz[0]);
            xyz[0] = xyz[0] + abc[0] * t.second; // координаты пересечения со сферой
            xyz[1] = xyz[1] + abc[1] * t.second;
            xyz[2] = xyz[2] + abc[2] * t.second;
            if (temp <= 0) {
                Cout_xyz(xyz); return 1;
            }
            if (curr_ht_ko == koef_osl[Lnum].begin()) {     // если летев внутрь частица сталкивается с поверхностью Земли                
                temp = Reflection(Lnum, xyz, abc, pp, type);
                if (temp) return -2;
            }
            else curr_ht_ko--; // частица летит вниз
        }
    }

    t.second = 0;
    if (curr_ht_ko == koef_osl[Lnum].end()) return -1;

    auto next = curr_ht_ko;
    next++;

    for (; next != koef_osl[Lnum].end(); ++curr_ht_ko) {
        R = next->first + 6371.0; // радиус сферы
        t = GetTequat(xyz, abc, R);
        if (t.first == 0) {
            return 0;
        }
        temp = l_opt - t.second * curr_ht_ko->second;
        if (temp <= 0) t.second = l_opt / curr_ht_ko->second;
        else l_opt = l_opt - t.second * curr_ht_ko->second;
        xyz[0] = xyz[0] + abc[0] * t.second; // координаты пересечения со сферой
        xyz[1] = xyz[1] + abc[1] * t.second;
        xyz[2] = xyz[2] + abc[2] * t.second;
        if (sqrt(xyz[2] * xyz[2] + xyz[1] * xyz[1] + xyz[0] * xyz[0]) - 6371 < -0.01)
            j = sqrt(xyz[2] * xyz[2] + xyz[1] * xyz[1] + xyz[0] * xyz[0]);
        if (temp <= 0) {
            Cout_xyz(xyz); return 1;
        }
        next++;
    }

    // если частица когда-либо оказалась в слое с нулевым коэффициентом, она вылетает за пределы атмосферы
    if (next == koef_osl[Lnum].end()) return -1;
}

/*bool ModelPerenosa::Reflection(double* xyz, double* abc) {
    // abc1 - вектор нормали к плоскости от которой отражается частица
    // abc2 = abc - направление до отражения
    // abc3 - направление после отражения

    // old basis - старый базис, основаный на abc1, вызывая функцию GetIzotr или GetLambert получаем координаты отраженного вектора abc3 в нем

    // если зеркальное отражение: abc3 = (0, sqrt(1-c3^2), c3)
    // c3 = cos teta = - a1*a2 - b1*b2 - c1*c2

    // new basis - новый базис
    // e1_old = (1,0,0)
    // e2_old = (0,1,0)
    // e3_old = (0,0,1)
    return 1;

    double abc2[3]{ abc[0], abc[1], abc[2] }; // задали abc2
    double abc3[3]{ 0, 0, 0 };
    double e1_old[3]{ 0, 0, 0 }; double e2_old[3]{ 0, 0, 0 }; double e3_old[3]{ 0, 0, 0 }; 
    double e1_new[3]{ 1, 0, 0 }; double e2_new[3]{ 0, 1, 0 }; double e3_new[3]{ 0, 0, 1 };  // задали новый базис

    double** a = new double* [3];   // матрица перехода к новому базису
    for (int i = 0; i < 3; i++)
        a[i] = new double[3];

    double v_length = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
    double abc1[3]{ xyz[0] / v_length, xyz[1] / v_length, xyz[2] / v_length }; // задали abc1

    // 1) находим координаты abc3 в старом базисе eсли отражение зеркальное
    bool zerkalnoe_otrazhenie = true;
    if (zerkalnoe_otrazhenie) {
        abc3[2] = -abc1[0] * abc2[0] - abc1[1] * abc2[1] - abc1[2] * abc2[2];
        abc3[1] = sqrt(1 - abc3[2]);
    }

    // если изотропное:
    // GetIzotr(abc);
    // abc3[0] = abc[0]; abc3[1] = abc[1]; abc3[2] = abc[2]; 
    // если Ламбертовское:
    // GetIzotr(abc);
    // abc3[0] = abc[0]; abc3[1] = abc[1]; abc3[2] = abc[2]; 

    // 2) находим e3_old
    e3_old[0] = abc1[0]; e3_old[1] = abc1[1]; e3_old[2] = abc1[2]; // задали e3 старого базиса

    // 3) находим e2_old. это пересечение двух плоскостей. первая - точка(0,0,0), abc1, abc2, вторая - перпендикулярна abc1
    // уравнение первой плоскости будет: A1 * x + B1 * y + C1 * z = 0
    // x abc1[0] abc2[0]
    // y abc1[1] abc2[1]
    // z abc1[2] abc2[2]
    double A1 = abc1[1] * abc2[2] - abc2[1] * abc1[2];
    double B1 = -abc1[0] * abc2[2] + abc2[0] * abc1[2];
    double C1 = abc1[0] * abc2[1] - abc2[0] * abc1[1];

    // уравнение второй плоскости будет: A2 * x + B2 * y + C2 * z = 0
    double A2 = abc1[0], B2 = abc1[1], C2 = abc1[2];

    // направляющий вектор прямой пересечения двух плоскостей:
    e2_old[0] = B1 * C2 - C1 * C2; e2_old[1] = C1 * A2 - A1 * C2; e2_old[2] = A1 * B2 - B1 * A2;

    // выбираем знак (чтобы вектор был по направлению abc2, а не -abc2)
    if ((e2_old[0] * abc2[0] + e2_old[1] * abc2[1] + e2_old[2] * abc2[2]) >= (-e2_old[0] * abc2[0] - e2_old[1] * abc2[1] - e2_old[2] * abc2[2])) {
        e2_old[0] = -e2_old[0]; e2_old[1] = -e2_old[1]; e2_old[2] = -e2_old[2];
    }

    // 4) находим e1_old
    // i         j          k
    // e2_old[0] e2_old[1] e2_old[2]
    // e3_old[0] e3_old[1] e3_old[2]
    e1_old[0] = e2_old[1] * e3_old[2] - e2_old[2] * e3_old[1];
    e1_old[1] = -e2_old[0] * e3_old[2] + e2_old[2] * e3_old[0];
    e1_old[2] = e2_old[0] * e3_old[1] - e2_old[1] * e3_old[0];

    findA(a[0], e1_old, e2_old, e3_old, e1_new);
    findA(a[1], e1_old, e2_old, e3_old, e2_new);
    findA(a[2], e1_old, e2_old, e3_old, e3_new);

    abc[0] = a[0][0] * e1_old[0] + a[0][1] * e2_old[0] + a[0][2] * e3_old[0];
    abc[1] = a[1][0] * e1_old[0] + a[1][1] * e2_old[0] + a[1][2] * e3_old[0];
    abc[2] = a[2][0] * e1_old[0] + a[2][1] * e2_old[0] + a[2][2] * e3_old[0];
}*/

// новая функция отражения

bool ModelPerenosa::Reflection(int Lnum, double* xyz, double* abc, double pp, int type) {
    // abc1 - вектор нормали к плоскости от которой отражается частица
    // abc2 = abc - направление до отражения
    // abc3 - направление после отражения
    if (EarthReflType(pp)) { 
        return 1; } // произошло поглощение
    
    double abc2[3]{ abc[0], abc[1], abc[2] }; // задали abc2

    if (type == 1) {
        double v_length = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
        double abc1[3]{ xyz[0] / v_length, xyz[1] / v_length, xyz[2] / v_length }; // задали abc1
        double spr = abc1[0] * abc2[0] + abc1[1] * abc2[1] + abc1[2] * abc2[2]; // скалярное произведение двух векторов - произведение их длин на косинус угла между ними, если есть координаты - сумма произведений соответствующих координат
        abc[0] = abc2[0] - 2 * spr * abc1[0];
        abc[1] = abc2[1] - 2 * spr * abc1[1];
        abc[2] = abc2[2] - 2 * spr * abc1[2];
        return 0; // произошло отражение
    }

    else {
        // old basis - старый базис, основаный на abc1, вызывая функцию GetIzotr или GetLambert получаем координаты отраженного вектора в нем
        // new basis - новый базис
        // e1_old = (1,0,0)
        // e2_old = (0,1,0)
        // e3_old = (0,0,1)

        double abc2[3]{ abc[0], abc[1], abc[2] }; // задали abc2
        double e1_old[3]{ 0, 0, 0 }; double e2_old[3]{ 0, 0, 0 }; double e3_old[3]{ 0, 0, 0 };
        double e1_new[3]{ 1, 0, 0 }; double e2_new[3]{ 0, 1, 0 }; double e3_new[3]{ 0, 0, 1 };  // задали новый базис - стандартный

        double** a = new double* [3];   // матрица перехода к новому базису
        for (int i = 0; i < 3; i++)
            a[i] = new double[3];

        double v_length = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
        double abc1[3]{ xyz[0] / v_length, xyz[1] / v_length, xyz[2] / v_length }; // задали abc1

        // если изотропное:
        if (type == 2) { GetIzotr(abc); }

        // если Ламбертовское:
        else { GetLambert(abc); }

        // 2) находим e3_old
        e3_old[0] = abc1[0]; e3_old[1] = abc1[1]; e3_old[2] = abc1[2]; // задали e3 старого базиса

        // 3) находим e2_old. это пересечение двух плоскостей. первая - точка(0,0,0), abc1, abc2, вторая - перпендикулярна abc1
        // уравнение первой плоскости будет: A1 * x + B1 * y + C1 * z = 0
        // x abc1[0] abc2[0]
        // y abc1[1] abc2[1]
        // z abc1[2] abc2[2]
        double A1 = abc1[1] * abc2[2] - abc2[1] * abc1[2];
        double B1 = -abc1[0] * abc2[2] + abc2[0] * abc1[2];
        double C1 = abc1[0] * abc2[1] - abc2[0] * abc1[1];

        // уравнение второй плоскости будет: A2 * x + B2 * y + C2 * z = 0
        double A2 = abc1[0], B2 = abc1[1], C2 = abc1[2];

        // направляющий вектор прямой пересечения двух плоскостей:
        e2_old[0] = B1 * C2 - C1 * C2; e2_old[1] = C1 * A2 - A1 * C2; e2_old[2] = A1 * B2 - B1 * A2;

        // выбираем знак (чтобы вектор был по направлению abc2, а не -abc2)
        if ((e2_old[0] * abc2[0] + e2_old[1] * abc2[1] + e2_old[2] * abc2[2]) >= (-e2_old[0] * abc2[0] - e2_old[1] * abc2[1] - e2_old[2] * abc2[2])) {
            e2_old[0] = -e2_old[0]; e2_old[1] = -e2_old[1]; e2_old[2] = -e2_old[2];
        }

        // 4) находим e1_old
        // i         j          k
        // e2_old[0] e2_old[1] e2_old[2]
        // e3_old[0] e3_old[1] e3_old[2]
        e1_old[0] = e2_old[1] * e3_old[2] - e2_old[2] * e3_old[1];
        e1_old[1] = -e2_old[0] * e3_old[2] + e2_old[2] * e3_old[0];
        e1_old[2] = e2_old[0] * e3_old[1] - e2_old[1] * e3_old[0];

        findA(a[0], e1_old, e2_old, e3_old, e1_new);
        findA(a[1], e1_old, e2_old, e3_old, e2_new);
        findA(a[2], e1_old, e2_old, e3_old, e3_new);

        abc[0] = a[0][0] * e1_old[0] + a[0][1] * e2_old[0] + a[0][2] * e3_old[0];
        abc[1] = a[1][0] * e1_old[0] + a[1][1] * e2_old[0] + a[1][2] * e3_old[0];
        abc[2] = a[2][0] * e1_old[0] + a[2][1] * e2_old[0] + a[2][2] * e3_old[0];
        
        for (int i = 0; i < 3; i++)
            delete[]a[i];
        delete[]a;
        return 0; // произошло отражение
    }
}
    
// решение методом Крамера
void ModelPerenosa::findA(double* a, double e1_old[], double e2_old[], double e3_old[], double e_new[])
{
    // вычисление определителя и проверка
    double det = e1_old[0] * e2_old[1] * e3_old[2] + e2_old[0] * e3_old[1] * e1_old[2] + e1_old[1] * e2_old[2] * e3_old[0] -
        e1_old[2] * e2_old[1] * e3_old[0] - e1_old[1] * e2_old[0] * e3_old[2] - e3_old[1] * e1_old[0] * e2_old[2];
    if (det == 0)
        return;

    double det1 = e_new[0] * e2_old[1] * e3_old[2] + e2_old[0] * e3_old[1] * e_new[2] + e_new[1] * e2_old[2] * e3_old[0] -
        e_new[2] * e2_old[1] * e3_old[0] - e_new[1] * e2_old[0] * e3_old[2] - e3_old[1] * e_new[0] * e2_old[2];

    double det2 = e1_old[0] * e_new[1] * e3_old[2] + e_new[0] * e3_old[1] * e1_old[2] + e1_old[1] * e_new[2] * e3_old[0] -
        e1_old[2] * e_new[1] * e3_old[0] - e1_old[1] * e_new[0] * e3_old[2] - e3_old[1] * e1_old[0] * e_new[2];

    double det3 = e1_old[0] * e2_old[1] * e_new[2] + e2_old[0] * e_new[1] * e1_old[2] + e1_old[1] * e2_old[2] * e_new[0] -
        e1_old[2] * e2_old[1] * e_new[0] - e1_old[1] * e2_old[0] * e_new[2] - e_new[1] * e1_old[0] * e2_old[2];

    a[0] = det1 / det; a[1] = det2 / det; a[2] = det3 / det;
}

std::pair<int,double> ModelPerenosa::GetTequat(double* xyz, double* abc, double R) {
    int xyz0[3]{ 0,0,0 };
    double D, x1, x2;
    double b = 2 * abc[0] * (xyz[0] - xyz0[0]) + 2 * abc[1] * (xyz[1] - xyz0[1]) + 2 * abc[2] * (xyz[2] - xyz0[2]);
    double a = abc[0] * abc[0] + abc[1] * abc[1] + abc[2] * abc[2];
    D = pow(b, 2) - 4 * a * (pow((xyz[0]-xyz0[0]),2) + pow((xyz[1] - xyz0[1]), 2) + pow((xyz[2] - xyz0[2]), 2) - R*R);
    if (D < 0) return std::make_pair(0, 0);
    else if (D == 0) {
        x1 = -b / (2 * a);
        if (x1 >= 0) return std::make_pair(1, x1); else return std::make_pair(0, 0);

    }
    else {
        x1 = (-b - sqrt(D)) / (2 * a);
        x2 = (-b + sqrt(D)) / (2 * a);
        if (x1 >= 0 && x2 >= 0) return std::make_pair(2, std::min(x1, x2));
        if (x1 >= 0 && x2 < 0) return std::make_pair(1, x1);
        if (x1 < 0 && x2 >= 0) return std::make_pair(1, x2);
        else return std::make_pair(0, 0);
    };
}

// проверка вылета из среды 
// вычисление координат очередной точки столкновения
/*double* ModelPerenosa::P3P4calcul(double* xyz, double* abc, double l, double* abc_prev) {
    if ((abc_prev[2] != abc[2]) & (abs(abc_prev[2]) != 0)) {
        double l1 = xyz[2] / abs(abc_prev[2]), l2 = l - l1;

        xyz[0] = xyz[0] + abc_prev[0] * l1;
        xyz[1] = xyz[1] + abc_prev[1] * l1;
        xyz[2] = xyz[2] + abc_prev[2] * l1;

        xyz[0] = xyz[0] + abc[0] * l2;
        xyz[1] = xyz[1] + abc[1] * l2;
        xyz[2] = xyz[2] + abc[2] * l2;
    }

    else {
        xyz[0] = xyz[0] + abc[0] * l;
        xyz[1] = xyz[1] + abc[1] * l;
        xyz[2] = xyz[2] + abc[2] * l;
    };

    return xyz;
}*/

// выбор типа столкновения (поглощение или рассеяние)
bool ModelPerenosa::P5type(int Lnum, std::vector<std::map<int, double>>& alb_rass, double* xyz) {
    double a = GetA();
    int curr_ht = xyz[2] - 6371, curr_ht_pos = -1;
    std::map<int, double>::iterator it = alb_rass[Lnum].begin();
    auto next = it;
    
    next++;
    for (;it != alb_rass[Lnum].end();) {
        if (curr_ht > next->first) {
            it++; next++;
        }
        else break;
    }

    if (a < it->second) {
        return 0;
    }// произошло рассеяние
    else return 1; // произошло поглощение
}

bool ModelPerenosa::EarthReflType(double pp) { // альбедо подстилающей поверхности
    double a = GetA();
    if (a <= pp) return 0; // произошло рассеяние
    else return 1; // произошло поглощение
}

// пересчет координат направления пробега
double* ModelPerenosa::P7napravl(float* angles, double** F, int Lnum, double* abc) {
    double* fi = new double[2], abct[3]; 
    double c = 1, m = 1;
    for (int i = 0; i < 3; i++)
        abct[i] = abc[i];
    while ((!((abc[2]>=-1)&((abc[2] <= 1)))) || (abs(c) == 1) || (abs(m)==1)) {
        m = getMa(angles, F, Lnum);
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

    /*std::cout << "x = " << xyz[0]
        << ", y = " << xyz[1]
        << ", z = " << xyz[2] << ";" << std::endl;*/
}

// функция для моделирования процесса переноса
int ModelPerenosa::ModPer(float* angles, double** F, int Lnum, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp, int type) {
    
    double* abc = new double[3], * xyz = new double[3];
    for (int i = 0; i < 3; i++)
        abc[i] = 0; 
    for (int i = 0; i < 2; i++)
        xyz[i] = 0;
    xyz[2] = 6371;
    int f;

    GetIzotr(abc);

    for (;;) {
        f = P2length(Lnum, koef_osl, alb_rass, xyz, abc, pp, type);

        if (f == -1)
        {
            // Произошел вылет за пределы среды через верхнюю границу
            Cout_xyz(xyz);
            CrossUp(abc[2]);
            delete[]abc;
            delete[]xyz;
            return 1;
        }
        //cout << "l = " << l << endl;

        if (f == -2)
        {
            // Произошло поглощение частицы поверхностью Земли
            Cout_xyz(xyz);
            CrossLow(abc[2]);
            delete[]abc;
            delete[]xyz;
            return 0;
        }

        if (P5type(Lnum, alb_rass, xyz)) {
            // Произошло поглощение
            Cout_xyz(xyz);
            delete[]abc;
            delete[]xyz;
            return 2;
        }

        abc = P7napravl(angles, F, Lnum, abc);

    }

}

int* ModelPerenosa::NModPer(int* t, float* angles, double** F, int Lnum, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp, int type)
{
    int k, j;
    for (int i = 0; i < 3; i++)
        t[i] = 0;
    for (int i = 0; i < kol; i++) {
        k = ModPer(angles, F, Lnum, koef_osl, alb_rass, pp, type);
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
        out << "Waves\tKtAbsorption\tKtUpCross\tUpCross\tKtLowCross\tLowCross" << std::endl;
        for (int i = 0; i < 5; i++) {
            out << waves[i] << '\t' << tBig[i][0] / (kol * 1.0) << '\t' << tBig[i][1] / (kol * 1.0) << '\t' << tBig[i][2];
            out << '\t' << tBig[i][3] / (kol * 1.0) << '\t' << tBig[i][4] << std::endl;
        }
    }
    out.close();
}

void ModelPerenosa::Modelirovanie(float* angles, double** F, double* waves, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp, int type)
{
    int* t = new int[3];
    double** tBig = new double* [5];
    for (int i = 0; i < 5; i++)
        tBig[i] = new double[5];
    for (int i = 0; i < 5; i++)
    {
        SetSum0();
        std::cout << "Данные для длины волны l=" << waves[i] << " мкм: " << std::endl << std::endl;
        t = NModPer(t, angles, F, i, koef_osl, alb_rass, pp, type);
        std::cout << "Произошло " << t[0] << " поглощений частиц поверхностью Земли. " << std::endl;
        std::cout << "Произошло " << t[1] << " вылетов за пределы среды через верхнюю границу. " << std::endl;
        std::cout << "Произошло " << t[2] << " поглощений. " << std::endl;

        std::cout << std::endl;
        CountK(t);

        std::cout << std::endl;
        std::cout << "Поток через верхнюю границу: " << GetSumUp() << std::endl;
        std::cout << "Поток через нижнюю границу: " << GetSumLow() << std::endl << std::endl << std::endl;
        
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
