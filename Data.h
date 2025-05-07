#pragma once

#include <iostream>
#include <fstream>
#include <time.h>
#include<vector>
#include <map>
#include "ModelPerenosa.h"

#define PI 3.14159265

class Data
{
public:
	Data() {};
	double* getWaves(double* waves);
	float* getM(float* angles);
	void getKoefOsl(std::vector<std::map<int, double>> &koef_osl, std::vector<std::map<int, double>> &alb_rass, double* waves);
	void GetMoleculScatterCoef(std::vector<std::map<int, double>>& mol_koef_scat, double* waves);
};

