#include "Trasser.h"

// Функции преобразования координат
void sphericalToCartesian(double r, double theta, double phi, double& x, double& y, double& z) 
{
    x = r * sin(theta) * cos(phi);
    y = r * sin(theta) * sin(phi);
    z = r * cos(theta);
}

double polar_angle(const double& x, const double& y)
{
    if (fabs(x) + fabs(y) < 0.000001)
    {
        return 0.0;
    }

    if (x < 0)
    {
        return atan(y / x) + 1.0 * M_PI;
    }
    else if (x > 0 && y >= 0)
    {
        return atan(y / x);
    }
    else if (x > 0 && y < 0)
    {
        return atan(y / x) + 2.0 * M_PI;
    }
    else if (y > 0 && x >= 0 && x <= 0)
    {
        return M_PI / 2.0;
    }
    else if (y < 0 && x >= 0 && x <= 0)
    {
        return  3.0 * M_PI / 2.0;
    }
    return 0.0;
}

void dekard_skorost(const double& z, const double& x, const double& y,
    const double& Vr, const double& Vphi, const double& Vtheta,
    double& Vz, double& Vx, double& Vy)
{
    double r_2, the_2, phi_2;
    r_2 = sqrt(x * x + y * y + z * z);
    the_2 = acos(z / r_2);
    phi_2 = polar_angle(x, y);

    Vx = Vr * sin(the_2) * cos(phi_2) + Vtheta * cos(the_2) * cos(phi_2) - Vphi * sin(phi_2);
    Vy = Vr * sin(the_2) * sin(phi_2) + Vtheta * cos(the_2) * sin(phi_2) + Vphi * cos(phi_2);
    Vz = Vr * cos(the_2) - Vtheta * sin(the_2);
}

void spherical_skorost(const double& z, const double& x, const double& y,
    const double& Vz, const double& Vx, const double& Vy,
    double& Vr, double& Vphi, double& Vtheta)
{
    double r_1, the_1, phi_1;

    r_1 = sqrt(x * x + y * y + z * z);
    the_1 = acos(z / r_1);
    phi_1 = polar_angle(x, y);

    Vr = Vx * sin(the_1) * cos(phi_1) + Vy * sin(the_1) * sin(phi_1) + Vz * cos(the_1);
    Vtheta = Vx * cos(the_1) * cos(phi_1) + Vy * cos(the_1) * sin(phi_1) - Vz * sin(the_1);
    Vphi = -Vx * sin(phi_1) + Vy * cos(phi_1);
}

void cartesianToSpherical(double x, double y, double z, double& r, double& theta, double& phi) 
{
    r = sqrt(x * x + y * y + z * z);
    theta = acos(z / r);
    phi = atan2(y, x);
    if (phi < 0) phi += 2 * M_PI;
}

// Функция трассировки одной магнитной линии (старая версия)
std::vector<std::vector<double>> traceFieldLine(
    double r_start, double theta_start, double phi_start,
    const std::vector<std::vector<Complex>>& B_lm,
    double R0, double Rss, int l_max,
    double dr, int max_steps) 
{

    std::vector<std::vector<double>> line;
    double r = r_start;
    double theta = theta_start;
    double phi = phi_start;

    // Добавляем начальную точку
    double x, y, z;
    sphericalToCartesian(r, theta, phi, x, y, z);
    double Br, Btheta, Bphi;
    double polar = 0.0;
    SphericalHarmonics::computeMagneticField(r, theta, phi, B_lm, R0, Rss, l_max, Br, Btheta, Bphi);
    if (Br > 0)
    {
        polar = 1.0;
    }
    else
    {
        polar = -1.0;
    }

    line.push_back({ x, y, z, polar });

    double k = 0.0;

    for (int step = 0; step < max_steps; step++) 
    {
        // Вычисляем магнитное поле в текущей точке
        SphericalHarmonics::computeMagneticField(r, theta, phi, B_lm, R0, Rss, l_max, Br, Btheta, Bphi);

        // Преобразуем компоненты поля в декартовы координаты
        double Bx = Br * sin(theta) * cos(phi) + Btheta * cos(theta) * cos(phi) - Bphi * sin(phi);
        double By = Br * sin(theta) * sin(phi) + Btheta * cos(theta) * sin(phi) + Bphi * cos(phi);
        double Bz = Br * cos(theta) - Btheta * sin(theta);

        // Нормируем вектор поля
        double B_mag = sqrt(Bx * Bx + By * By + Bz * Bz);
        if (B_mag < 1e-12) break;

        Bx /= B_mag;
        By /= B_mag;
        Bz /= B_mag;

        // Делаем шаг вдоль поля
        double ds = dr * polar; // Шаг по длине дуги
        x += ds * Bx;
        y += ds * By;
        z += ds * Bz;

        // Преобразуем обратно в сферические координаты
        cartesianToSpherical(x, y, z, r, theta, phi);

        // Проверяем условия остановки
        if (r > Rss || r < R0) break;

        k = -1.0;
        if (Br > 0) k = 1.0;

        // Добавляем точку к линии
        line.push_back({ x, y, z, k });
    }

    return line;
}

// Новая оптимизированная функция трассировки с использованием предвычисленной сетки
std::vector<std::vector<double>> traceFieldLineWithGrid(
    double r_start, double theta_start, double phi_start,
    const MagneticFieldGrid& grid,
    double dr, int max_steps) 
{
    std::vector<std::vector<double>> line;
    double r = r_start;
    double theta = theta_start;
    double phi = phi_start;

    // Добавляем начальную точку
    double x, y, z;
    sphericalToCartesian(r, theta, phi, x, y, z);
    
    // Получаем магнитное поле из сетки
    auto field = grid.interpolateField(r, theta, phi);
    double Br = field.Br;
    double Btheta = field.Btheta;
    double Bphi = field.Bphi;
    
    double polar = (Br > 0) ? 1.0 : -1.0;
    line.push_back({ x, y, z, polar });

    for (int step = 0; step < max_steps; step++) 
    {
        // Получаем магнитное поле из интерполированной сетки
        auto field = grid.interpolateField(r, theta, phi);
        Br = field.Br;
        Btheta = field.Btheta;
        Bphi = field.Bphi;

        // Преобразуем компоненты поля в декартовы координаты
        double Bx = Br * sin(theta) * cos(phi) + Btheta * cos(theta) * cos(phi) - Bphi * sin(phi);
        double By = Br * sin(theta) * sin(phi) + Btheta * cos(theta) * sin(phi) + Bphi * cos(phi);
        double Bz = Br * cos(theta) - Btheta * sin(theta);

        // Нормируем вектор поля
        double B_mag = sqrt(Bx * Bx + By * By + Bz * Bz);
        if (B_mag < 1e-12) break;

        Bx /= B_mag;
        By /= B_mag;
        Bz /= B_mag;

        // Делаем шаг вдоль поля
        double ds = dr * polar; // Шаг по длине дуги
        x += ds * Bx;
        y += ds * By;
        z += ds * Bz;

        // Преобразуем обратно в сферические координаты
        cartesianToSpherical(x, y, z, r, theta, phi);

        // Проверяем условия остановки
        if (r > grid.getRss() || r < grid.getR0()) break;

        double k = (Br > 0) ? 1.0 : -1.0;
        line.push_back({ x, y, z, k });
    }

    return line;
}

// Основная функция для генерации всех магнитных линий
void generateMagneticFieldLines(
    const std::vector<std::vector<Complex>>& B_lm,
    const std::vector<std::vector<double>>& PHI,
    const std::vector<std::vector<double>>& THETA,
    const std::vector<std::vector<double>>& Br_2d,
    double R0, double Rss, int l_max,
    const std::string& filename) 
{

    std::ofstream file("open_" + filename);
    std::ofstream file2("closed_" + filename);
    file << "VARIABLES = X, Y, Z, I" << std::endl;
    file2 << "VARIABLES = X, Y, Z, I" << std::endl;

    std::ofstream file3("coronal_phole_" + filename);
    file3 << "VARIABLES = phi, the, I" << std::endl;

    double dr = 0.005 * R0; // Шаг трассировки  0.005
    int max_steps = 100000;  // Максимальное количество шагов

    int line_count = 0;
    int ik = 0;
    short int step_i = 1;
    // Выбираем стартовые точки на фотосфере (каждую 10-ю точку по сетке)
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < PHI.size(); i += step_i)
    {
        #pragma omp critical (ergertet4) 
        {
            ik++;
            if (ik % 10 == 0)
            {
                cout << "Step " << ik << "    from  " << PHI.size() / step_i << endl;
            }
        }
        for (int j = 0; j < PHI[0].size(); j += 1) 
        {
            double theta = THETA[i][j];
            double phi = PHI[i][j];
            double Br = Br_2d[i][j];

            auto line = traceFieldLine(R0, theta, phi, B_lm, R0, Rss, l_max, dr, max_steps);

            double x1, y1, z1;
            x1 = line[line.size() - 1][0];
            y1 = line[line.size() - 1][1];
            z1 = line[line.size() - 1][2];

            double r = sqrt(x1 * x1 + y1 * y1 + z1 * z1);

            if (r > Rss * 0.99)
            {
                // Записываем линию в файл
                #pragma omp critical (ge) 
                {
                    file << "ZONE T=\"Line" << line_count++ << "\" I=" << line.size() << " F=POINT" << std::endl;
                    for (const auto& point : line)
                    {
                        file << point[0] << " " << point[1] << " " << point[2] << " " << point[3] << std::endl;
                    }

                    file3 << phi << " " << theta << " " << 1.0 << endl;
                }
            }
            else
            {
                // Записываем линию в файл
                #pragma omp critical (ge) 
                {
                    file2 << "ZONE T=\"Line" << line_count++ << "\" I=" << line.size() << " F=POINT" << std::endl;
                    for (const auto& point : line)
                    {
                        file2 << point[0] << " " << point[1] << " " << point[2] << " " << point[3] << std::endl;
                    }

                    file3 << phi << " " << theta << " " << 0.0 << endl;
                }
            }
        }
    }

    file.close();
    file2.close();
    file3.close();
    std::cout << "Generated " << line_count << " magnetic field lines in " << filename << std::endl;
}


void generate_Coronal_hole(
    const std::vector<std::vector<Complex>>& B_lm,
    const std::vector<std::vector<double>>& PHI,
    const std::vector<std::vector<double>>& THETA,
    const std::vector<std::vector<double>>& Br_2d,
    double R0, double Rss, int l_max,
    const std::string& filename)
{

    std::ofstream file(filename);
    file << "VARIABLES = phi, the, I" << std::endl;

    double dr = 0.005 * R0; // Шаг трассировки
    int max_steps = 10000;  // Максимальное количество шагов

    int line_count = 0;

    // Выбираем стартовые точки на фотосфере (каждую 10-ю точку по сетке)
    for (int i = 0; i < PHI.size(); i += 1)
    {
        for (int j = 0; j < PHI[0].size(); j += 1)
        {
            double theta = THETA[i][j];
            double phi = PHI[i][j];
            double Br = Br_2d[i][j];

            auto line = traceFieldLine(R0, theta, phi, B_lm, R0, Rss, l_max, dr, max_steps);

            double x1, y1, z1;
            x1 = line[line.size() - 1][0];
            y1 = line[line.size() - 1][1];
            z1 = line[line.size() - 1][2];

            double r = sqrt(x1 * x1 + y1 * y1 + z1 * z1);

            if (r > Rss * 0.99)
            {
                file << phi << " " << theta << " " << 1.0 << endl;
            }
            else
            {
                file << phi << " " << theta << " " << 0.0 << endl;
            }
        }
    }

    file.close();
}

// Новая оптимизированная функция для корональных дыр с предвычисленной сеткой
void generate_Coronal_hole_optimized(
    const std::vector<std::vector<double>>& PHI,
    const std::vector<std::vector<double>>& THETA,
    const std::vector<std::vector<double>>& Br_2d,
    const MagneticFieldGrid& grid,
    const std::string& filename)
{
    std::ofstream file(filename);
    file << "VARIABLES = phi, the, I" << std::endl;

    double dr = 0.005 * grid.getR0(); // Шаг трассировки
    int max_steps = 10000;  // Максимальное количество шагов

    std::cout << "Starting optimized coronal hole generation..." << std::endl;
    std::cout << "Processing " << PHI.size() << " x " << PHI[0].size() << " = " 
              << (PHI.size() * PHI[0].size()) << " points" << std::endl;

    int total_points = PHI.size() * PHI[0].size();
    int processed_points = 0;
    int progress_step = total_points / 20; // Показываем прогресс каждые 5%

    // Выбираем стартовые точки на фотосфере
    for (int i = 0; i < PHI.size(); i += 1)
    {
        for (int j = 0; j < PHI[0].size(); j += 1)
        {
            double theta = THETA[i][j];
            double phi = PHI[i][j];
            double Br = Br_2d[i][j];

            // Используем новую оптимизированную функцию трассировки
            auto line = traceFieldLineWithGrid(grid.getR0(), theta, phi, grid, dr, max_steps);

            double x1, y1, z1;
            x1 = line[line.size() - 1][0];
            y1 = line[line.size() - 1][1];
            z1 = line[line.size() - 1][2];

            double r = sqrt(x1 * x1 + y1 * y1 + z1 * z1);

            if (r > grid.getRss() * 0.99)
            {
                file << phi << " " << theta << " " << 1.0 << endl;
            }
            else
            {
                file << phi << " " << theta << " " << 0.0 << endl;
            }

            processed_points++;
            if (progress_step > 0 && processed_points % progress_step == 0) {
                int progress = (100 * processed_points) / total_points;
                std::cout << "Tracing progress: " << progress << "%" << std::endl;
            }
        }
    }

    file.close();
    std::cout << "Optimized coronal hole generation completed!" << std::endl;
}

void generateSphere(double radius, int numTheta, int numPhi,
    std::vector<Point3D>& vertices,
    std::vector<std::vector<int>>& connectivity) {

    vertices.clear();
    connectivity.clear();

    // Генерация вершин
    for (int i = 0; i <= numTheta; ++i) {
        double theta = static_cast<double>(i) / numTheta * M_PI; // от 0 до ?

        for (int j = 0; j <= numPhi; ++j) {
            double phi = static_cast<double>(j) / numPhi * 2.0 * M_PI; // от 0 до 2?

            Point3D p;
            p.x = radius * sin(theta) * cos(phi);
            p.y = radius * sin(theta) * sin(phi);
            p.z = radius * cos(theta);

            vertices.push_back(p);
        }
    }

    // Генерация коннективити (четырехугольники)
    for (int i = 0; i < numTheta; ++i) {
        for (int j = 0; j < numPhi; ++j) {
            int current = i * (numPhi + 1) + j;
            int next = (i + 1) * (numPhi + 1) + j;

            std::vector<int> quad(4);
            quad[0] = current + 1;       // P1
            quad[1] = current;           // P2
            quad[2] = next;              // P3
            quad[3] = next + 1;          // P4

            connectivity.push_back(quad);
        }
    }
}

void writeTecplotFile(const std::string& filename,
    const std::vector<Point3D>& vertices,
    const std::vector<std::vector<int>>& connectivity) {

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка открытия файла!" << std::endl;
        return;
    }

    // Заголовок
    file << "TITLE = Sphere\n";
    file << "VARIABLES = X, Y, Z\n";
    file << "ZONE T=\"Sphere\", N=" << vertices.size()
        << ", E=" << connectivity.size()
        << ", F=FEPOINT, ET=QUADRILATERAL\n";

    // Вершины
    file << std::fixed << std::setprecision(6);
    for (const auto& vertex : vertices) {
        file << std::setw(12) << vertex.x
            << std::setw(12) << vertex.y
            << std::setw(12) << vertex.z << "\n";
    }

    // Коннективити (индексы начинаются с 1)
    for (const auto& quad : connectivity) {
        file << std::setw(8) << quad[0] + 1
            << std::setw(8) << quad[1] + 1
            << std::setw(8) << quad[2] + 1
            << std::setw(8) << quad[3] + 1 << "\n";
    }

    file.close();
}