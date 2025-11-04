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
#include <omp.h>
#include <chrono>
#include <thread>
#include <mutex>

struct Point3D
{
    double x, y, z;
};

// ‘орвардные декларации всех классов и структур
class PFSSData;
class SphericalHarmonics;
class MagneticFieldGrid;
struct SphericalHarmonicData;
struct MagneticFieldPoint;



using Complex = std::complex<double>;


//#define M_PI 3.14159265


using namespace std;

#include "PFSS_Data.h"
#include "Trasser.h"
#include "SphericalHarmonics.h"
#include "MagneticFieldGrid.h"