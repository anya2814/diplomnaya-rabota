#pragma once

#include <iostream>
#include <fstream>
#include <time.h>
#include<vector>
#include <cmath>
#include <string>
#include <map>

#define PI 3.14159265359

struct KDistData {
	double lambda_um = -1.0;                 // средн€€ длина волны (мкм)
	double Ci[4]{ 0,0,0,0 };                   // C1..C4
	std::vector<std::map<int, double>> k_layers; // [m][z_key] = k_m(z), m=0..3
};

class Data
{
public:
	Data() {};
	double* getWaves(double* waves);
	float* getM(float* angles);
	void getKoefOsl(std::vector<std::map<int, double>> &koef_osl, std::vector<std::map<int, double>> &alb_rass, double* waves);
	void getMoleculScatterCoef(std::vector<std::map<int, double>>& mol_koef_scat, double* waves);
	void getMoleculConsCoef(std::vector<std::map<int, double>>& mol_cons, double* waves);
	bool readMolAbsKdistFile(const std::string& filename, KDistData& out);
};

