#include"Header.h"



// Структура для хранения сферических гармоник и их производных
struct SphericalHarmonicData {
    Complex Y;
    Complex dY_dtheta;
    Complex dY_dphi;
};


// Функция для вычисления сферических гармоник и их производных
// Старафай функция не работает при нулевых углах
void computeSphericalHarmonicsWithDerivatives_old(
    int l_max,
    double theta,
    double phi,
    std::vector<std::vector<SphericalHarmonicData>>& Y_data) {

    double x = cos(theta);
    double sin_theta = sin(theta);

    size_t array_size = gsl_sf_legendre_array_n(l_max);
    std::vector<double> legendre_array(array_size);
    std::vector<double> legendre_deriv_array(array_size);

    gsl_sf_legendre_deriv_array(
        GSL_SF_LEGENDRE_SPHARM,
        l_max,
        x,
        legendre_array.data(),
        legendre_deriv_array.data()
    );

    // ВЫЧИСЛЯЕМ С ФАЗОВЫМ МНОЖИТЕЛЕМ (-1)^m
    for (int l = 0; l <= l_max; l++) {
        for (int m = -l; m <= l; m++) {
            int abs_m = std::abs(m);
            int idx = gsl_sf_legendre_array_index(l, abs_m);

            double legendre_val = legendre_array[idx];
            double legendre_deriv = legendre_deriv_array[idx];

            // Фазовый множитель Condon-Shortley (-1)^m
            double phase = (abs_m % 2 == 0) ? 1.0 : -1.0;

            Complex Y, dY_dtheta, dY_dphi;

            if (m >= 0) {
                Y = phase * legendre_val * std::exp(Complex(0, m * phi));
                dY_dtheta = phase * (-sin_theta * legendre_deriv) * std::exp(Complex(0, m * phi));
                dY_dphi = phase * legendre_val * Complex(0, m) * std::exp(Complex(0, m * phi));
            }
            else {
                // Для отрицательных m используем свойство симметрии
                Complex Y_pos = phase * legendre_val * std::exp(Complex(0, abs_m * phi));
                Complex dY_dtheta_pos = phase * (-sin_theta * legendre_deriv) * std::exp(Complex(0, abs_m * phi));
                Complex dY_dphi_pos = phase * legendre_val * Complex(0, abs_m) * std::exp(Complex(0, abs_m * phi));

                Y = std::pow(-1.0, abs_m) * std::conj(Y_pos);
                dY_dtheta = std::pow(-1.0, abs_m) * std::conj(dY_dtheta_pos);
                dY_dphi = std::pow(-1.0, abs_m) * std::conj(dY_dphi_pos);
            }

            Y_data[l][m + l] = { Y, dY_dtheta, dY_dphi };
        }
    }
}

void computeSphericalHarmonicsWithDerivatives(
    int l_max,
    double theta,
    double phi,
    std::vector<std::vector<SphericalHarmonicData>>& Y_data) {

    double x = cos(theta);
    double sin_theta = sin(theta);

    // Обработка особых случаев на полюсах
    if (std::abs(x - 1.0) < 1e-12 || std::abs(x + 1.0) < 1e-12) {
        // На полюсах все сферические гармоники с m ? 0 равны 0
        // А производные на полюсах требуют специальной обработки

        for (int l = 0; l <= l_max; l++) {
            for (int m = -l; m <= l; m++) {
                int m_index = m + l;

                if (m == 0) {
                    // Только для m = 0 сферические гармоники ненулевые на полюсах
                    double legendre_val = gsl_sf_legendre_sphPlm(l, 0, x);
                    Complex Y = legendre_val; // ?-зависимость отсутствует при m=0

                    // На полюсах производные по ? вычисляются аналитически
                    Complex dY_dtheta, dY_dphi;

                    if (l == 0) {
                        dY_dtheta = 0.0;
                    }
                    else {
                        // Для l > 0 и m=0: dY_l0/d? = -sqrt((l(l+1))/(4?)) * Y_l1
                        // Но на полюсах Y_l1 = 0, поэтому нам нужно использовать предел
                        if (std::abs(x - 1.0) < 1e-12) { // Северный полюс
                            dY_dtheta = 0.0;
                        }
                        else { // Южный полюс
                            dY_dtheta = 0.0;
                        }
                    }

                    dY_dphi = 0.0; // При m=0 производная по ? всегда 0

                    Y_data[l][m_index] = { Y, dY_dtheta, dY_dphi };
                }
                else {
                    // Для m ? 0 на полюсах все компоненты равны 0
                    Y_data[l][m_index] = { 0.0, 0.0, 0.0 };
                }
            }
        }
        return;
    }

    // Стандартное вычисление для неособых точек
    size_t array_size = gsl_sf_legendre_array_n(l_max);
    std::vector<double> legendre_array(array_size);
    std::vector<double> legendre_deriv_array(array_size);

    gsl_sf_legendre_deriv_array(
        GSL_SF_LEGENDRE_SPHARM,
        l_max,
        x,
        legendre_array.data(),
        legendre_deriv_array.data()
    );

    // Вычисляем с фазовым множителем (-1)^m
    for (int l = 0; l <= l_max; l++) {
        for (int m = -l; m <= l; m++) {
            int abs_m = std::abs(m);
            int idx = gsl_sf_legendre_array_index(l, abs_m);

            double legendre_val = legendre_array[idx];
            double legendre_deriv = legendre_deriv_array[idx];

            // Фазовый множитель Condon-Shortley (-1)^m
            double phase = (abs_m % 2 == 0) ? 1.0 : -1.0;

            Complex Y, dY_dtheta, dY_dphi;

            if (m >= 0) {
                Y = phase * legendre_val * std::exp(Complex(0, m * phi));
                dY_dtheta = phase * (-sin_theta * legendre_deriv) * std::exp(Complex(0, m * phi));
                dY_dphi = phase * legendre_val * Complex(0, m) * std::exp(Complex(0, m * phi));
            }
            else {
                // Для отрицательных m используем свойство симметрии
                Complex Y_pos = phase * legendre_val * std::exp(Complex(0, abs_m * phi));
                Complex dY_dtheta_pos = phase * (-sin_theta * legendre_deriv) * std::exp(Complex(0, abs_m * phi));
                Complex dY_dphi_pos = phase * legendre_val * Complex(0, abs_m) * std::exp(Complex(0, abs_m * phi));

                Y = std::pow(-1.0, abs_m) * std::conj(Y_pos);
                dY_dtheta = std::pow(-1.0, abs_m) * std::conj(dY_dtheta_pos);
                dY_dphi = std::pow(-1.0, abs_m) * std::conj(dY_dphi_pos);
            }

            Y_data[l][m + l] = { Y, dY_dtheta, dY_dphi };
        }
    }
}

// Функция для вычисления магнитного поля в любой точке
void computeMagneticField(
    double r, double theta, double phi,
    const std::vector<std::vector<Complex>>& a_lm,
    double R0, double Rss, int l_max,
    double& Br, double& Btheta, double& Bphi) {

    // Вычисляем сферические гармоники и их производные
    std::vector<std::vector<SphericalHarmonicData>> Y_data(l_max + 1);
    for (int l = 0; l <= l_max; l++) {
        Y_data[l].resize(2 * l + 1);
    }

    computeSphericalHarmonicsWithDerivatives(l_max, theta, phi, Y_data);

    Br = 0.0;
    Btheta = 0.0;
    Bphi = 0.0;

    // ОТЛАДОЧНЫЙ ВЫВОД
    //std::cout << "Debug: r=" << r << ", R0=" << R0 << ", Rss=" << Rss << std::endl;

    for (int l = 0; l <= l_max; l++) {
        for (int m = -l; m <= l; m++) {
            int m_index = m + l;
            Complex a = a_lm[l][m_index];

            const auto& Y = Y_data[l][m_index];

            // ПРАВИЛЬНАЯ формула для радиальной компоненты
            double term1 = (l) * std::pow(r / R0, l - 1);
            double term2 = (l + 1) * std::pow(r / R0, -(l + 2)) * std::pow(Rss/R0, 2 * l + 1);
            double denominator = l + (l + 1) * std::pow(Rss/R0, 2 * l + 1);

            double radial_Br = (term1 + term2) / denominator;

            // Радиальная функция для поперечных компонентов
            double radial_Btrans = (std::pow(R0 / r, l + 2) - std::pow(r / R0, l - 1) * std::pow(R0 / Rss, 2 * l + 1)) / denominator;

            // Компоненты магнитного поля
            Complex br_term = a * radial_Br * Y.Y;
            Br += std::real(br_term);

            Complex btheta_term = a * radial_Btrans * Y.dY_dtheta;
            Btheta += std::real(btheta_term);

            // Аккуратная обработка полюсов
            if (std::abs(sin(theta)) > 1e-12) {
                Complex bphi_term = a * radial_Btrans * (1.0 / sin(theta)) * Y.dY_dphi;
                Bphi += std::real(bphi_term);
            }

            if (false)
            {
                // ОТЛАДОЧНЫЙ ВЫВОД для первых нескольких гармоник
                if (l <= 2 && std::abs(a) > 1e-10) {
                    std::cout << "l=" << l << ", m=" << m
                        << ": radial_Br=" << radial_Br
                        << " (term1=" << term1 << ", term2=" << term2
                        << ", den=" << denominator << ")"
                        << ", a_lm=" << a
                        << ", Y=" << Y.Y
                        << ", product=" << br_term
                        << std::endl;
                }
            }
        }
    }

    // Дополнительная проверка: сравнение с прямой реконструкцией
    double Br_direct = 0.0;
    for (int l = 0; l <= l_max; l++) {
        for (int m = -l; m <= l; m++) {
            int m_index = m + l;
            Complex a = a_lm[l][m_index];
            const auto& Y = Y_data[l][m_index];
            Br_direct += std::real(a * Y.Y);
        }
    }

    /*std::cout << "Br (PFSS) = " << Br << std::endl;
    std::cout << "Br (direct) = " << Br_direct << std::endl;
    std::cout << "Btheta = " << Btheta << std::endl;
    std::cout << "Bphi = " << Bphi << std::endl;*/
}


// Разложение Br на фотосфере - нахождение коэффициентов разложения
std::vector<std::vector<Complex>> computeSphericalHarmonics_photosphere(
    const std::vector<std::vector<double>>& Br_2d,
    const std::vector<std::vector<double>>& PHI,
    const std::vector<std::vector<double>>& THETA,
    int l_max = 15) {

    int n_lon = Br_2d.size();
    int n_lat = Br_2d[0].size();

    std::vector<std::vector<Complex>> B_lm(l_max + 1);
    for (int l = 0; l <= l_max; l++) {
        B_lm[l].resize(2 * l + 1, 0.0);
    }

    double dphi = 2.0 * M_PI / n_lon;
    double dtheta = M_PI / (n_lat - 1);  // Исправлено: обычно n_lat включает полюса

    for (int i = 0; i < n_lon; i++) {
        for (int j = 0; j < n_lat; j++) {
            double phi = PHI[i][j];
            double theta = THETA[i][j];
            double sin_theta = sin(theta);

            // Пропускаем полюса, где sin_theta = 0
            if (sin_theta < 1e-12) continue;

            double Br = Br_2d[i][j];

            for (int l = 0; l <= l_max; l++) {
                for (int m = 0; m <= l; m++) {
                    double legendre = gsl_sf_legendre_sphPlm(l, m, cos(theta));
                    Complex Y_lm = legendre * std::exp(Complex(0, m * phi));

                    // Умножаем на sin(theta) - элемент площади сферы
                    B_lm[l][m + l] += Br * std::conj(Y_lm) * sin_theta * dphi * dtheta;
                }
            }
        }
    }

    // Заполняем отрицательные m с правильной индексацией
    for (int l = 0; l <= l_max; l++) {
        for (int m = 1; m <= l; m++) {
            // Индекс для отрицательного m: l - m
            B_lm[l][l - m] = std::pow(-1.0, m) * std::conj(B_lm[l][l + m]);
        }
    }

    return B_lm;
}


double reconstructBr(
    double theta,
    double phi,
    const std::vector<std::vector<Complex>>& B_lm,
    int l_max = -1) {

    // Если l_max не указан, используем максимальный доступный
    if (l_max == -1) {
        l_max = B_lm.size() - 1;
    }

    // Вычисляем сферические гармоники и их производные
    /*std::vector<std::vector<SphericalHarmonicData>> Y_data(l_max + 1);
    for (int l = 0; l <= l_max; l++) {
        Y_data[l].resize(2 * l + 1);
    }
    computeSphericalHarmonicsWithDerivatives(l_max, theta, phi, Y_data);*/

    double Br = 0.0;
    double Br2 = 0.0;

    for (int l = 0; l <= l_max; l++) {
        // Проверяем, что для данного l есть коэффициенты
        if (l >= B_lm.size()) continue;

        for (int m_index = 0; m_index < B_lm[l].size(); m_index++) {
            int m = m_index - l;  // Преобразуем индекс в значение m

            // Вычисляем сферическую гармонику Y_lm(?, ?)
            double legendre = gsl_sf_legendre_sphPlm(l, std::abs(m), cos(theta));
            Complex Y_lm;

            if (m >= 0) {
                Y_lm = legendre * std::exp(Complex(0, m * phi));
            }
            else {
                // Для отрицательных m используем свойство симметрии
                Y_lm = std::pow(-1.0, std::abs(m)) * legendre * std::exp(Complex(0, m * phi));
            }

            

            // Добавляем вклад этого коэффициента
            Br += std::real(B_lm[l][m_index] * Y_lm);

            //auto Y = Y_data[l][m_index];
            //Br2 += std::real(B_lm[l][m_index] * Y.Y);
        }
    }

    //cout << "Br " << Br << "  " << Br2 << endl;

    return Br;
}

void writeComparisonToFile(
    const std::string& filename,
    const std::vector<std::vector<double>>& PHI,
    const std::vector<std::vector<double>>& THETA,
    const std::vector<std::vector<double>>& Br_orig,
    const std::vector<std::vector<Complex>>& B_lm, const double& RR,
    int step_lon = 1,
    int step_lat = 1, int l_max = 15) {

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Записываем заголовок
    file << "TITLE = HP  VARIABLES = phi, the, Br_orig, Br_my, Br, Bthe, Brphi" << std::endl;

    // Записываем данные
    file << std::scientific << std::setprecision(10);

    for (int i = 0; i < PHI.size(); i += step_lon) {
        for (int j = 0; j < PHI[0].size(); j += step_lat) {
            double phi = PHI[i][j];
            double theta = THETA[i][j];
            double Br_original = Br_orig[i][j];
            double Br_reconstructed = reconstructBr(theta, phi, B_lm);
            double Br, Btheta, Bphi;

            //cout << "1  " << theta << " " << phi << endl;
            computeMagneticField(RR, theta, phi, B_lm, 1.0, 2.5, l_max, Br, Btheta, Bphi);
            //cout << "2   " << Br  << endl;

            file << phi << " " << theta << " "
                << Br_original << " " << Br_reconstructed << " " << Br << " " << Btheta << " " << Bphi << std::endl;
        }
    }

    file.close();
    std::cout << "Comparison data written to " << filename << std::endl;
}

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

        auto B_lm = computeSphericalHarmonics_photosphere(data.Br_2d, data.PHI_2d, data.THETA_2d, l_max);

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
                    double Br_reconstructed = reconstructBr(theta, phi, B_lm);

                    std::cout << "The=" << theta << ", Phi=" << phi
                        << ": original=" << Br_original
                        << ", reconstruct=" << Br_reconstructed
                        << ", error=" << std::abs(Br_original - Br_reconstructed)
                        << std::endl;
                }
            }
        }

        cout << "To file" << endl;
        //writeComparisonToFile("1_comparison_50.txt", data.PHI_2d, data.THETA_2d, data.Br_2d, B_lm, 1.0, 1, 1, l_max);
        //writeComparisonToFile("2.5_comparison_50.txt", data.PHI_2d, data.THETA_2d, data.Br_2d, B_lm, 2.5, 1, 1, l_max);

        // Вычисляем поле в конкретной точке
        /*double r = 1.0, theta = M_PI / 3, phi = M_PI / 8;
        double Br, Btheta, Bphi;

        computeMagneticField(r, theta, phi, B_lm, 1.0, 2.5, l_max, Br, Btheta, Bphi);
        std::cout << "B = (" << Br << ", " << Btheta << ", " << Bphi << ")" << std::endl;
        cout << reconstructBr(theta, phi, B_lm) << endl;*/

        generateMagneticFieldLines(B_lm, data.PHI_2d, data.THETA_2d, data.Br_2d, 1.0, 2.5, l_max, "magnetic_lines.txt");
        generate_Coronal_hole(B_lm, data.PHI_2d, data.THETA_2d, data.Br_2d, 1.0, 2.5, l_max, "coronal_phole.txt");


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