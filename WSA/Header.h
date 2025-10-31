#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <complex>

using namespace std;
using Complex = std::complex<double>;

struct Point3D {
    double x, y, z;
};


#define M_PI 3.14159265

class PFSSData;

#include "PFSS_Data.h"
#include "Trasser.h"


void computeMagneticField(
    double r, double theta, double phi,
    const std::vector<std::vector<Complex>>& a_lm,
    double R0, double Rss, int l_max,
    double& Br, double& Btheta, double& Bphi);