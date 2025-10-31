#pragma once

#include"Header.h"

void sphericalToCartesian(double r, double theta, double phi, double& x, double& y, double& z);
void cartesianToSpherical(double x, double y, double z, double& r, double& theta, double& phi);

std::vector<std::vector<double>> traceFieldLine(
    double r_start, double theta_start, double phi_start,
    const std::vector<std::vector<Complex>>& B_lm,
    double R0, double Rss, int l_max,
    double dr, int max_steps);

void generateMagneticFieldLines(
    const std::vector<std::vector<Complex>>& B_lm,
    const std::vector<std::vector<double>>& PHI,
    const std::vector<std::vector<double>>& THETA,
    const std::vector<std::vector<double>>& Br_2d,
    double R0, double Rss, int l_max,
    const std::string& filename);

void writeTecplotFile(const std::string& filename,
    const std::vector<Point3D>& vertices,
    const std::vector<std::vector<int>>& connectivity);

void generateSphere(double radius, int numTheta, int numPhi,
    std::vector<Point3D>& vertices,
    std::vector<std::vector<int>>& connectivity);