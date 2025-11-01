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

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using Complex = std::complex<double>;

struct Point3D {
    double x, y, z;
};

#define M_PI 3.14159265

// ‘орвардные декларации всех классов и структур
class PFSSData;
class SphericalHarmonics;
class MagneticFieldGrid;
struct SphericalHarmonicData;
struct MagneticFieldPoint;

#include "PFSS_Data.h"
#include "Trasser.h"
#include "SphericalHarmonics.h"
#include "MagneticFieldGrid.h"