#pragma once

#include "Header.h"

// Структура для хранения сферических гармоник и их производных
struct SphericalHarmonicData {
    Complex Y;
    Complex dY_dtheta;
    Complex dY_dphi;
};

class SphericalHarmonics {
public:
    // Функция для вычисления сферических гармоник и их производных
    // Старая функция не работает при нулевых углах
    static void computeSphericalHarmonicsWithDerivatives_old(
        int l_max,
        double theta,
        double phi,
        std::vector<std::vector<SphericalHarmonicData>>& Y_data);

    // Улучшенная функция для вычисления сферических гармоник и их производных
    static void computeSphericalHarmonicsWithDerivatives(
        int l_max,
        double theta,
        double phi,
        std::vector<std::vector<SphericalHarmonicData>>& Y_data);

    // Разложение Br на фотосфере - нахождение коэффициентов разложения
    static std::vector<std::vector<Complex>> computeSphericalHarmonics_photosphere(
        const std::vector<std::vector<double>>& Br_2d,
        const std::vector<std::vector<double>>& PHI,
        const std::vector<std::vector<double>>& THETA,
        int l_max = 15);

    // Восстановление Br по коэффициентам сферических гармоник
    static double reconstructBr(
        double theta,
        double phi,
        const std::vector<std::vector<Complex>>& B_lm,
        int l_max = -1);

    // Функция для вычисления магнитного поля в любой точке
    static void computeMagneticField(
        double r, double theta, double phi,
        const std::vector<std::vector<Complex>>& a_lm,
        double R0, double Rss, int l_max,
        double& Br, double& Btheta, double& Bphi);

    // Запись сравнения в файл
    static void writeComparisonToFile(
        const std::string& filename,
        const std::vector<std::vector<double>>& PHI,
        const std::vector<std::vector<double>>& THETA,
        const std::vector<std::vector<double>>& Br_orig,
        const std::vector<std::vector<Complex>>& B_lm,
        const double& RR,
        int step_lon = 1,
        int step_lat = 1,
        int l_max = 15);
};