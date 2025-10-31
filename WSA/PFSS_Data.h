#pragma once


#include"Header.h"


class PFSSData {
public:
    int n_lon, n_lat;
    std::vector<std::vector<double>> Br_2d;
    std::vector<std::vector<double>> PHI_2d;
    std::vector<std::vector<double>> THETA_2d;
    std::vector<double> phi_1d;
    std::vector<double> theta_1d;

    // Методы для чтения данных
    bool loadMetadata(const std::string& filename);
    bool load2DArray(const std::string& filename, std::vector<std::vector<double>>& array, int n_rows, int n_cols);
    bool load1DArray(const std::string& filename, std::vector<double>& array);

    // Основной метод загрузки всех данных
    bool loadAllData(const std::string& base_filename);

    // Вспомогательные методы
    void printInfo() const;
    double getBr(int i, int j) const { return Br_2d[i][j]; }
    double getPhi(int i, int j) const { return PHI_2d[i][j]; }
    double getTheta(int i, int j) const { return THETA_2d[i][j]; }
};
