#pragma once

#include <iostream>
#include <fstream>
#include <time.h>
#include<vector>
#include "ModelPerenosa.h"

#define PI 3.14159265

class Data
{
public:
	Data() {};
	double* getWaves(double* waves);
	double** getF(double** F, float* mass);
	float* getM(float* mass);
	double** getd(double** d, double* waves);
};

