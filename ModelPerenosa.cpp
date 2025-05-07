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

double ModelPerenosa::GetWeight(double* xyz, double* abc) {
    double res = (xyz[0] * abc[0] + xyz[1] * abc[1] + xyz[2] * abc[2]) / sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
    return abs(res);
}

// учет пересечений верхней площадки с весом 1/|(ns, w)|
void ModelPerenosa::CrossUp(double* xyz, double* abc)
{
    sumUp = sumUp + 1.0 / (kol * GetWeight(xyz, abc));
}

double ModelPerenosa::GetSumUp()
{
    return sumUp;
}

// учет пересечений нижней площадки с весом 1/|(ns, w)|
void ModelPerenosa::CrossLow(double* xyz, double* abc)
{
    sumLow = sumLow + 1.0 / (kol * GetWeight(xyz, abc));
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

// НОВАЯ ФУНКЦИЯ выбор длины свободного пробега l + проверка вылета из среды 
// вычисление координат очередной точки столкновения
int ModelPerenosa::P2length(int Lnum, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double* xyz, double* abc, double pp, int type) {
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


// функция отражения
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
        double spr = abc1[0] * abc2[0] + abc1[1] * abc2[1] + abc1[2] * abc2[2]; // скалярное произведение двух векторов - произведение их длин на косинус угла между ними, 
                                                                                // если есть координаты - сумма произведений соответствующих координат
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
        double e1_old[3]{ 1, 0, 0 }; double e2_old[3]{ 0, 1, 0 }; double e3_old[3]{ 0, 0, 1 };
        double e1_new[3]{ 0, 0, 0 }; double e2_new[3]{ 0, 0, 0 }; double e3_new[3]{ 0, 0, 0 };  // задали новый базис - стандартный

        double *old_coord = new double[3];

        double** a = new double* [3];   // матрица перехода к новому базису
        for (int i = 0; i < 3; i++)
            a[i] = new double[3];

        double v_length = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
        double abc1[3]{ xyz[0] / v_length, xyz[1] / v_length, xyz[2] / v_length }; // задали abc1

        // если изотропное:
        if (type == 2) {
            GetIzotr(old_coord); 
        }

        // если Ламбертовское:
        else { GetLambert(old_coord); }

        // 2) находим e3_new
        e3_new[0] = abc1[0]; e3_new[1] = abc1[1]; e3_new[2] = abc1[2]; // задали e3 нового базиса

        // 3) находим e2_new. это пересечение двух плоскостей. первая - точка(0,0,6371), abc1, abc2, вторая - перпендикулярна abc1
        // уравнение первой плоскости будет: A1 * x + B1 * y + C1 * z = 0
        
        // проверка abc2 на параллельность вектору abc1
        while (abs(abc1[0] * abc2[0] + abc1[1] * abc2[1] + abc1[2] * abc2[2]) > 0.99) {
            GetIzotr(abc2); // генерируем вектор равномерно распределенный в одной полуплоскости
            // с вероятностью 1/2 умножаем вектор на минус 1 чтобы был равномерно распределен в двух полуплоскостях, то есть в сфере
            if (GetA() < 0.5) { abc2[0] = -abc2[0]; abc2[1] = -abc2[1]; abc2[2] = -abc2[2]; }
        }
        
        // x abc1[0] abc2[0]
        // y abc1[1] abc2[1]
        // z abc1[2] abc2[2]
        double A1 = abc1[1] * abc2[2] - abc2[1] * abc1[2];
        double B1 = -abc1[0] * abc2[2] + abc2[0] * abc1[2];
        double C1 = abc1[0] * abc2[1] - abc2[0] * abc1[1];

        // уравнение второй плоскости будет: A2 * x + B2 * y + C2 * z = 0
        double A2 = abc1[0], B2 = abc1[1], C2 = abc1[2];

        // направляющий вектор прямой пересечения двух плоскостей:
        e2_new[0] = B1 * C2 - C1 * C2; e2_new[1] = C1 * A2 - A1 * C2; e2_new[2] = A1 * B2 - B1 * A2;

        // выбираем знак (чтобы вектор был по направлению abc2, а не -abc2)
        if ((e2_new[0] * abc2[0] + e2_new[1] * abc2[1] + e2_new[2] * abc2[2]) < (-e2_new[0] * abc2[0] - e2_new[1] * abc2[1] - e2_new[2] * abc2[2])) {
            e2_new[0] = -e2_new[0]; e2_new[1] = -e2_new[1]; e2_new[2] = -e2_new[2];
        }
        
        abc2[0] = abc[0]; abc2[1] = abc[1]; abc2[2] = abc[2]; // на случай если брали случайный вектор для построения плоскости

        // нормировка
        v_length = sqrt(e2_new[0] * e2_new[0] + e2_new[1] * e2_new[1] + e2_new[2] * e2_new[2]);
        for (int i = 0; i < 3; i++)
        {
            e2_new[i] = e2_new[i] / v_length;
        }

        // 4) находим e1_new
        // i         j          k
        // e2_new[0] e2_new[1] e2_new[2]
        // e3_new[0] e3_new[1] e3_new[2]
        e1_new[0] = e2_new[1] * e3_new[2] - e2_new[2] * e3_new[1];
        e1_new[1] = -e2_new[0] * e3_new[2] + e2_new[2] * e3_new[0];
        e1_new[2] = e2_new[0] * e3_new[1] - e2_new[1] * e3_new[0];

        // e1_new = a[0][0] * e1_old + a[1][0] * e2_old + a[2][0] * e3_old
        // e2_new = a[0][1] * e1_old + a[1][1] * e2_old + a[2][1] * e3_old
        // e3_new = a[0][2] * e1_old + a[1][2] * e2_old + a[2][2] * e3_old

        findA(a[0], e1_old, e2_old, e3_old, e1_new);
        findA(a[1], e1_old, e2_old, e3_old, e2_new);
        findA(a[2], e1_old, e2_old, e3_old, e3_new);

        abc[0] = a[0][0] * old_coord[0] + a[1][0] * old_coord[1] + a[2][0] * old_coord[2];
        abc[1] = a[0][1] * old_coord[0] + a[1][1] * old_coord[1] + a[2][1] * old_coord[2];
        abc[2] = a[0][2] * old_coord[0] + a[1][2] * old_coord[1] + a[2][2] * old_coord[2];
        int h = 0;

        for (int i = 0; i < 3; i++)
            delete[]a[i];
        delete[]a;
        delete[]old_coord;
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
int ModelPerenosa::ModPer(float* angles, double** F, int Lnum, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp, int type) {
    
    double* abc = new double[3], * xyz = new double[3];
    for (int i = 0; i < 3; i++)
        abc[i] = 0; 
    for (int i = 0; i < 2; i++)
        xyz[i] = 0;
    xyz[2] = 6371;
    int f;

    GetIzotr(abc);

    for (;;) {
        f = P2length(Lnum, mol_koef_rass, koef_osl, alb_rass, xyz, abc, pp, type);

        if (f == -1)
        {
            // Произошел вылет за пределы среды через верхнюю границу
            CrossUp(xyz, abc);
            delete[]abc;
            delete[]xyz;
            return 1;
        }
        //cout << "l = " << l << endl;

        if (f == -2)
        {
            // Произошло поглощение частицы поверхностью Земли
            CrossLow(xyz, abc);
            delete[]abc;
            delete[]xyz;
            return 0;
        }

        if (P5type(Lnum, alb_rass, xyz)) {
            // Произошло поглощение
            delete[]abc;
            delete[]xyz;
            return 2;
        }

        abc = P7napravl(angles, F, Lnum, abc);

    }

}

int* ModelPerenosa::NModPer(int* t, float* angles, double** F, int Lnum, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp, int type)
{
    int k, j;
    for (int i = 0; i < 3; i++)
        t[i] = 0;
    for (int i = 0; i < kol; i++) {
        k = ModPer(angles, F, Lnum, mol_koef_rass, koef_osl, alb_rass, pp, type);
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
void ModelPerenosa::OutToFile(double** tBig, double* waves, double pp)
{
    for (int i = 0; i < 5; i++) {
        std::ofstream out;        
        out.open("Results_wave_" + std::to_string(waves[i]).substr(0, 5) + ".txt", std::ios::app);      // открываем файл для записи
        if (out.is_open())
        {
            out << pp << '\t' << tBig[i][0] / (kol * 1.0) << '\t' << tBig[i][1] / (kol * 1.0) << '\t' << tBig[i][2];
            out << '\t' << tBig[i][3] / (kol * 1.0) << '\t' << tBig[i][4] << std::endl;
        };
        out.close();
    }
}

void ModelPerenosa::Modelirovanie(float* angles, double* waves, std::vector<std::map<int, double>>& mol_koef_rass, std::vector<std::map<int, double>>& koef_osl, std::vector<std::map<int, double>>& alb_rass, double pp, int type)
{
    int* t = new int[3];
    double** tBig = new double* [5];
    for (int i = 0; i < 5; i++)
        tBig[i] = new double[5];
    double** F = new double* [N+1];       // массив значений эмпирической функции распределения направлений рассеяния
    for (int i = 0; i < N+1; i++)
        F[i] = new double[23];

    for (int i = 0; i < 5; i++)
    {
        SetSum0();
        F = getF(F, angles, i, mol_koef_rass[i], koef_osl[i], alb_rass[i]);      // функция распределения угла рассеяния
        std::cout << "Данные для длины волны l=" << waves[i] << " мкм: " << std::endl << std::endl;
        t = NModPer(t, angles, F, i, mol_koef_rass, koef_osl, alb_rass, pp, type);
        /*std::cout << "Произошло " << t[0] << " поглощений частиц поверхностью Земли. " << std::endl;
        std::cout << "Произошло " << t[1] << " вылетов за пределы среды через верхнюю границу. " << std::endl;
        std::cout << "Произошло " << t[2] << " поглощений. " << std::endl;*/

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

    OutToFile(tBig, waves, pp);

    for (int i = 0; i < N + 1; i++)
        delete[]F[i];
    delete[]F;
    delete[]t;
    for (int i = 0; i < 5; i++)
        delete[]tBig[i];
    delete[]tBig;
}

// вычисление F по формуле трапеций
double** ModelPerenosa::getF(double** F, float* mass, int Lnum, std::map<int, double>& mol_koef_rass, std::map<int, double>& koef_osl, std::map<int, double>& alb_rass)
{
    double* ind = new double [N];
    double sum = 0;
    double one = 1;

    int i = 0, col;
    if (F[0][0] == -1) {
        for (std::map<int, double>::iterator it = mol_koef_rass.begin(); it != mol_koef_rass.end(); it++) {
            F[0][i] = it->first;
            i++;
        }
    }

    std::ifstream H;
    float read;
    
    std::map<int, double>::iterator it_mol = mol_koef_rass.begin();
    std::map<int, double>::iterator next_mol = it_mol; next_mol++;
    std::map<int, double>::iterator it_aer = koef_osl.begin();
    std::map<int, double>::iterator next_aer = it_aer; next_aer++; 
    std::map<int, double>::iterator it_aer_scat = alb_rass.begin();
    std::map<int, double>::iterator next_aer_scat = it_aer_scat; next_aer_scat++;

    H.open("Indikatrisa.txt");
    if (H.is_open()) {
        for (int i = 0; i < 5; i++) {
            H >> read;
        }
        bool flag = false;

        for (int i = 0; i < N; i++)
        {
            H >> read;

            col = 0;
            for (int j = 0; j < 5; j++) {
                H >> read;

                if (j == Lnum) {
                    ind[i] = read / (2 * PI);
                    for (; next_mol != mol_koef_rass.end();) {
                        ind[i] = ind[i] * (it_aer->second) * (it_aer_scat->second) + 
                            3 / 8 * (1 + one * mass[i] * mass[i]) * (it_mol->second) / ((it_aer->second) * (it_aer_scat->second) + (it_mol->second));
                        if (flag)
                            sum = sum + (ind[i] + ind[i - 1]) / 2.0 * abs(mass[i] - mass[i - 1]);
                        F[i+1][col] = sum;
                        col++;
                        if (next_aer->first == next_mol->first) {
                            it_aer++;
                            next_aer++;
                        }
                        if (next_aer_scat->first == next_mol->first) {
                            it_aer_scat++;
                            next_aer_scat++;
                        }
                        it_mol++;
                        next_mol++;
                    }
                    it_mol = mol_koef_rass.begin(); next_mol = it_mol; next_mol++;
                    it_aer = koef_osl.begin(); next_aer = it_aer; next_aer++;
                    it_aer_scat = alb_rass.begin(); next_aer_scat = it_aer_scat; next_aer_scat++;
                }
            }
            flag = true;
        }
        H.close();
    }

    delete[]ind;

    return F;
}
