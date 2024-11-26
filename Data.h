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
	double** getF(double** F, float* mass);
	float* getM(float* mass);
	void getKoefOsl(std::vector<std::map<int, double>> &koef_osl, std::vector<std::map<int, double>> &alb_rass, double* waves);
};

