#include"Header.h"

// Основная программа
int main() {
    PFSSData data;

    std::cout << "Loading solar magnetic field data..." << std::endl;

    if (data.loadAllData("solar_data")) {
        std::cout << "Data loaded successfully!" << std::endl;
        data.printInfo();

        short int l_max = 50;

        // Здесь можно добавить вызовы для PFSS модели
        // Например: data.calculatePFSS();

        auto B_lm = SphericalHarmonics::computeSphericalHarmonics_photosphere(data.Br_2d, data.PHI_2d, data.THETA_2d, l_max);

        // Вывводим коэффициенты
        if (false)
        {
            std::cout << "Prin coefficient:  " << B_lm.size() << std::endl;

            std::cout << "Coefficients B_lm:" << std::endl;
            std::cout << "Number of l values: " << B_lm.size() << std::endl;

            for (int l = 0; l < B_lm.size(); l++)
            {
                std::cout << "l = " << l << " (m from -" << l << " to " << l << "):" << std::endl;
                for (int m_index = 0; m_index < B_lm[l].size(); m_index++)
                {
                    int m = m_index - l;  // Преобразуем индекс в значение m
                    std::cout << "  B_" << l << "," << m << " = " << B_lm[l][m_index]
                        << " (magnitude: " << std::abs(B_lm[l][m_index]) << ")" << std::endl;
                }
                std::cout << std::endl;
            }
        }

        if (false)
        {
            // Проверка восстановления в точках сетки
            std::cout << "Proverka:" << std::endl;
            for (int i = 0; i < data.PHI_2d.size(); i += 10) {  // Каждую 10-ю точку
                for (int j = 0; j < data.THETA_2d[0].size(); j += 10) {
                    double theta = data.THETA_2d[i][j];
                    double phi = data.PHI_2d[i][j];
                    double Br_original = data.Br_2d[i][j];
                    double Br_reconstructed = SphericalHarmonics::reconstructBr(theta, phi, B_lm);

                    std::cout << "The=" << theta << ", Phi=" << phi
                        << ": original=" << Br_original
                        << ", reconstruct=" << Br_reconstructed
                        << ", error=" << std::abs(Br_original - Br_reconstructed)
                        << std::endl;
                }
            }
        }

        cout << "To file" << endl;
        //SphericalHarmonics::writeComparisonToFile("1_comparison_50.txt", data.PHI_2d, data.THETA_2d, data.Br_2d, B_lm, 1.0, 1, 1, l_max);
        //SphericalHarmonics::writeComparisonToFile("2.5_comparison_50.txt", data.PHI_2d, data.THETA_2d, data.Br_2d, B_lm, 2.5, 1, 1, l_max);

        // Вычисляем поле в конкретной точке
        /*double r = 1.0, theta = M_PI / 3, phi = M_PI / 8;
        double Br, Btheta, Bphi;

        SphericalHarmonics::computeMagneticField(r, theta, phi, B_lm, 1.0, 2.5, l_max, Br, Btheta, Bphi);
        std::cout << "B = (" << Br << ", " << Btheta << ", " << Bphi << ")" << std::endl;
        cout << SphericalHarmonics::reconstructBr(theta, phi, B_lm) << endl;*/

        //generateMagneticFieldLines(B_lm, data.PHI_2d, data.THETA_2d, data.Br_2d, 1.0, 2.5, l_max, "magnetic_lines.txt");
        
        // НОВАЯ ОПТИМИЗИРОВАННАЯ ВЕРСИЯ с предвычисленной сеткой
        std::cout << "\nCreating optimized magnetic field grid..." << std::endl;
        
        // Параметры сетки: разрешение по углам как в исходных данных, настраиваемое разрешение по радиусу
        int nr_grid = 50;  // Количество слоев по радиусу (настраиваемый параметр)
        int ntheta_grid = data.THETA_2d[0].size();  // Берем разрешение по theta из исходных данных
        int nphi_grid = data.PHI_2d.size();         // Берем разрешение по phi из исходных данных
        
        MagneticFieldGrid grid(1.0, 2.5, nr_grid, ntheta_grid, nphi_grid);
        grid.initializeGrid(B_lm, l_max);
        
        std::cout << "Grid created successfully! Starting optimized coronal hole calculation..." << std::endl;
        generate_Coronal_hole_optimized(data.PHI_2d, data.THETA_2d, data.Br_2d, grid, "coronal_phole_optimized.txt");

        if (false)
        {
            // Генерируем сферу в техплоте
            // Параметры сферы
            double radius = 1.0;      // единичный радиус
            int numTheta = 200;        // количество разбиений по theta (меридианы)
            int numPhi = 200;          // количество разбиений по phi (параллели)

            std::vector<Point3D> vertices;
            std::vector<std::vector<int>> connectivity;

            // Генерация сферы
            generateSphere(radius, numTheta, numPhi, vertices, connectivity);

            // Запись в файл
            writeTecplotFile("sphere.txt", vertices, connectivity);
        }

    }
    else {
        std::cerr << "Failed to load data!" << std::endl;
        return 1;
    }

    return 0;
}