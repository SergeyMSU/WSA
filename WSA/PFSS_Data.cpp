#include "PFSS_Data.h"

bool PFSSData::loadMetadata(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open metadata file " << filename << std::endl;
        return false;
    }

    std::string key;
    double value;

    while (file >> key >> value) {
        if (key == "n_lon") n_lon = static_cast<int>(value);
        else if (key == "n_lat") n_lat = static_cast<int>(value);
        // Остальные метаданные пока не используем, но можем прочитать если нужно
    }

    file.close();
    std::cout << "Metadata loaded: n_lon=" << n_lon << ", n_lat=" << n_lat << std::endl;
    return true;
}

bool PFSSData::load2DArray(const std::string& filename, std::vector<std::vector<double>>& array, int n_rows, int n_cols) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }

    // Инициализируем массив
    array.resize(n_rows, std::vector<double>(n_cols));

    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_cols; ++j) {
            if (!(file >> array[i][j])) {
                std::cerr << "Error reading data from " << filename << " at position (" << i << "," << j << ")" << std::endl;
                return false;
            }
        }
    }

    file.close();
    std::cout << "2D array loaded from " << filename << " with shape (" << n_rows << "," << n_cols << ")" << std::endl;
    return true;
}

bool PFSSData::load1DArray(const std::string& filename, std::vector<double>& array) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }

    double value;
    while (file >> value) {
        array.push_back(value);
    }

    file.close();
    std::cout << "1D array loaded from " << filename << " with size " << array.size() << std::endl;
    return true;
}

bool PFSSData::loadAllData(const std::string& base_filename) {
    // Загружаем метаданные
    if (!loadMetadata(base_filename + "_metadata.txt")) {
        return false;
    }

    // Загружаем 2D массивы
    if (!load2DArray(base_filename + "_Br_2d.dat", Br_2d, n_lon, n_lat)) return false;
    if (!load2DArray(base_filename + "_PHI_2d.dat", PHI_2d, n_lon, n_lat)) return false;
    if (!load2DArray(base_filename + "_THETA_2d.dat", THETA_2d, n_lon, n_lat)) return false;

    // Загружаем 1D массивы
    if (!load1DArray(base_filename + "_phi_1d.dat", phi_1d)) return false;
    if (!load1DArray(base_filename + "_theta_1d.dat", theta_1d)) return false;

    return true;
}

void PFSSData::printInfo() const {
    std::cout << "\n=== PFSS Data Info ===" << std::endl;
    std::cout << "Dimensions: " << n_lon << " x " << n_lat << std::endl;
    std::cout << "Br_2d range: " << Br_2d[0][0] << " to " << Br_2d[n_lon - 1][n_lat - 1] << std::endl;
    std::cout << "PHI range: " << PHI_2d[0][0] << " to " << PHI_2d[n_lon - 1][n_lat - 1] << " radians" << std::endl;
    std::cout << "THETA range: " << THETA_2d[0][0] << " to " << THETA_2d[n_lon - 1][n_lat - 1] << " radians" << std::endl;
    std::cout << "phi_1d size: " << phi_1d.size() << std::endl;
    std::cout << "theta_1d size: " << theta_1d.size() << std::endl;

    // Проверим несколько точек
    std::cout << "\nSample points:" << std::endl;
    std::cout << "Br_2d[100][100] = " << getBr(100, 100) << std::endl;
    std::cout << "PHI[100][100] = " << getPhi(100, 100) << " rad = " << getPhi(100, 100) * 180 / M_PI << " deg" << std::endl;
    std::cout << "THETA[100][100] = " << getTheta(100, 100) << " rad = " << getTheta(100, 100) * 180 / M_PI << " deg" << std::endl;
}



