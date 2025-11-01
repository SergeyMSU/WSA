#include "MagneticFieldGrid.h"

// Конструктор
MagneticFieldGrid::MagneticFieldGrid(double R0, double Rss, int nr, int ntheta, int nphi)
    : R0(R0), Rss(Rss), nr(nr), ntheta(ntheta), nphi(nphi) {
    
    // Инициализация координатных сеток
    r_grid.resize(nr);
    theta_grid.resize(ntheta);
    phi_grid.resize(nphi);
    
    // Создание логарitmической сетки по радиусу
    double log_R0 = std::log(R0);
    double log_Rss = std::log(Rss);
    for (int i = 0; i < nr; i++) {
        double log_r = log_R0 + (log_Rss - log_R0) * i / (nr - 1);
        r_grid[i] = std::exp(log_r);
    }
    
    // Равномерная сетка по theta (от 0 до PI)
    for (int i = 0; i < ntheta; i++) {
        theta_grid[i] = M_PI * i / (ntheta - 1);
    }
    
    // Равномерная сетка по phi (от 0 до 2*PI)
    for (int i = 0; i < nphi; i++) {
        phi_grid[i] = 2.0 * M_PI * i / nphi; // Периодическая по phi
    }
    
    // Инициализация сетки магнитного поля
    field_grid.resize(nr);
    for (int i = 0; i < nr; i++) {
        field_grid[i].resize(ntheta);
        for (int j = 0; j < ntheta; j++) {
            field_grid[i][j].resize(nphi);
        }
    }
}

// Инициализация сетки с предвычислением магнитного поля (с OpenMP)
void MagneticFieldGrid::initializeGrid(const std::vector<std::vector<Complex>>& B_lm, int l_max) {
    std::cout << "Initializing magnetic field grid..." << std::endl;
    std::cout << "Grid size: " << nr << " x " << ntheta << " x " << nphi 
              << " = " << (nr * ntheta * nphi) << " points" << std::endl;
    
#ifdef _OPENMP
    std::cout << "Using OpenMP with " << omp_get_max_threads() << " threads" << std::endl;
#endif
    
    int total_points = nr * ntheta * nphi;
    int computed_points = 0;
    int progress_step = total_points / 20;
    
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) reduction(+:computed_points)
#endif
    for (int ir = 0; ir < nr; ir++) {
        int local_count = 0;
        for (int it = 0; it < ntheta; it++) {
            for (int ip = 0; ip < nphi; ip++) {
                double r = r_grid[ir];
                double theta = theta_grid[it];
                double phi = phi_grid[ip];
                
                double Br, Btheta, Bphi;
                SphericalHarmonics::computeMagneticField(r, theta, phi, B_lm, R0, Rss, l_max, 
                                                       Br, Btheta, Bphi);
                
                field_grid[ir][it][ip] = MagneticFieldPoint(Br, Btheta, Bphi);
                local_count++;
            }
        }
        computed_points += local_count;
        
#ifdef _OPENMP
        if (omp_get_thread_num() == 0 && computed_points % progress_step == 0) {
#else
        if (computed_points % progress_step == 0) {
#endif
            int progress = (100 * computed_points) / total_points;
            std::cout << "Progress: " << progress << "%" << std::endl;
        }
    }
    
    std::cout << "Grid initialization completed!" << std::endl;
}

// Вспомогательная функция для поиска индексов и весов интерполяции
void MagneticFieldGrid::findGridIndices(double r, double theta, double phi, 
                                       int& ir, int& itheta, int& iphi,
                                       double& wr, double& wtheta, double& wphi) const {
    
    // Поиск индекса по r (логарithмический поиск)
    double log_r = std::log(r);
    double log_R0 = std::log(R0);
    double log_Rss = std::log(Rss);
    double r_normalized = (log_r - log_R0) / (log_Rss - log_R0);
    double r_index_real = r_normalized * (nr - 1);
    
    ir = static_cast<int>(std::floor(r_index_real));
    ir = std::max(0, std::min(ir, nr - 2));
    wr = r_index_real - ir;
    
    // Поиск индекса по theta
    double theta_index_real = theta / M_PI * (ntheta - 1);
    itheta = static_cast<int>(std::floor(theta_index_real));
    itheta = std::max(0, std::min(itheta, ntheta - 2));
    wtheta = theta_index_real - itheta;
    
    // Поиск индекса по phi (с учетом периодичности)
    double phi_normalized = phi;
    if (phi_normalized < 0) phi_normalized += 2.0 * M_PI;
    if (phi_normalized >= 2.0 * M_PI) phi_normalized -= 2.0 * M_PI;
    
    double phi_index_real = phi_normalized / (2.0 * M_PI) * nphi;
    iphi = static_cast<int>(std::floor(phi_index_real));
    iphi = iphi % nphi; // Обеспечиваем периодичность
    wphi = phi_index_real - std::floor(phi_index_real);
}

// Интерполяция магнитного поля в произвольной точке
MagneticFieldPoint MagneticFieldGrid::interpolateField(double r, double theta, double phi) const {
    if (!isInGrid(r, theta, phi)) {
        return MagneticFieldPoint(0.0, 0.0, 0.0);
    }
    
    int ir, itheta, iphi;
    double wr, wtheta, wphi;
    findGridIndices(r, theta, phi, ir, itheta, iphi, wr, wtheta, wphi);
    
    // Трилинейная интерполяция
    // Получаем 8 соседних точек
    auto get_field = [&](int dr, int dt, int dp) -> MagneticFieldPoint {
        int ir_idx = std::min(ir + dr, nr - 1);
        int it_idx = std::min(itheta + dt, ntheta - 1);
        int ip_idx = (iphi + dp) % nphi; // Периодичность по phi
        return field_grid[ir_idx][it_idx][ip_idx];
    };
    
    // Интерполяция по r
    auto interp_r = [&](int dt, int dp) -> MagneticFieldPoint {
        auto f0 = get_field(0, dt, dp);
        auto f1 = get_field(1, dt, dp);
        return MagneticFieldPoint(
            f0.Br * (1 - wr) + f1.Br * wr,
            f0.Btheta * (1 - wr) + f1.Btheta * wr,
            f0.Bphi * (1 - wr) + f1.Bphi * wr
        );
    };
    
    // Интерполяция по theta
    auto interp_rt = [&](int dp) -> MagneticFieldPoint {
        auto f0 = interp_r(0, dp);
        auto f1 = interp_r(1, dp);
        return MagneticFieldPoint(
            f0.Br * (1 - wtheta) + f1.Br * wtheta,
            f0.Btheta * (1 - wtheta) + f1.Btheta * wtheta,
            f0.Bphi * (1 - wtheta) + f1.Bphi * wtheta
        );
    };
    
    // Интерполяция по phi
    auto f0 = interp_rt(0);
    auto f1 = interp_rt(1);
    
    return MagneticFieldPoint(
        f0.Br * (1 - wphi) + f1.Br * wphi,
        f0.Btheta * (1 - wphi) + f1.Btheta * wphi,
        f0.Bphi * (1 - wphi) + f1.Bphi * wphi
    );
}

// Проверка, находится ли точка в пределах сетки
bool MagneticFieldGrid::isInGrid(double r, double theta, double phi) const {
    return (r >= R0 && r <= Rss && theta >= 0 && theta <= M_PI);
}