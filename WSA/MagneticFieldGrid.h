#pragma once

#include "Header.h"

// Структура для хранения компонентов магнитного поля в точке
struct MagneticFieldPoint {
    double Br, Btheta, Bphi;
    
    MagneticFieldPoint() : Br(0.0), Btheta(0.0), Bphi(0.0) {}
    MagneticFieldPoint(double br, double bt, double bp) : Br(br), Btheta(bt), Bphi(bp) {}
};

// Класс для работы с предвычисленной сеткой магнитного поля
class MagneticFieldGrid {
private:
    // Параметры сетки
    double R0, Rss;
    int nr, ntheta, nphi;
    
    // Координаты сетки
    std::vector<double> r_grid;
    std::vector<double> theta_grid;
    std::vector<double> phi_grid;
    
    // Предвычисленные значения магнитного поля
    // Индексация: [ir][itheta][iphi]
    std::vector<std::vector<std::vector<MagneticFieldPoint>>> field_grid;
    
    // Вспомогательные функции для поиска индексов
    void findGridIndices(double r, double theta, double phi, 
                        int& ir, int& itheta, int& iphi,
                        double& wr, double& wtheta, double& wphi) const;

public:
    // Конструктор
    MagneticFieldGrid(double R0, double Rss, int nr, int ntheta, int nphi);
    
    // Инициализация сетки с предвычислением магнитного поля
    void initializeGrid(const std::vector<std::vector<Complex>>& B_lm, int l_max);
    
    // Интерполяция магнитного поля в произвольной точке
    MagneticFieldPoint interpolateField(double r, double theta, double phi) const;
    
    // Получение параметров сетки
    double getR0() const { return R0; }
    double getRss() const { return Rss; }
    int getNr() const { return nr; }
    int getNtheta() const { return ntheta; }
    int getNphi() const { return nphi; }
    
    // Проверка, находится ли точка в пределах сетки
    bool isInGrid(double r, double theta, double phi) const;
};