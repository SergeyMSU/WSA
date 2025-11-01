#include "Trasser.h"

// Функции преобразования координат
void sphericalToCartesian(double r, double theta, double phi, double& x, double& y, double& z) 
{
    x = r * sin(theta) * cos(phi);
    y = r * sin(theta) * sin(phi);
    z = r * cos(theta);
}

void cartesianToSpherical(double x, double y, double z, double& r, double& theta, double& phi) 
{
    r = sqrt(x * x + y * y + z * z);
    theta = acos(z / r);
    phi = atan2(y, x);
    if (phi < 0) phi += 2 * M_PI;
}


// Функция трассировки одной магнитной линии
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

    double dr = 0.005 * R0; // Шаг трассировки
    int max_steps = 10000;  // Максимальное количество шагов

    int line_count = 0;

    // Выбираем стартовые точки на фотосфере (каждую 10-ю точку по сетке)
    for (int i = 0; i < PHI.size(); i += 10) 
    {
        for (int j = 0; j < PHI[0].size(); j += 10) 
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
                file << "ZONE T=\"Line" << line_count++ << "\" I=" << line.size() << " F=POINT" << std::endl;
                for (const auto& point : line)
                {
                    file << point[0] << " " << point[1] << " " << point[2] << " " << point[3] << std::endl;
                }
            }
            else
            {
                // Записываем линию в файл
                file2 << "ZONE T=\"Line" << line_count++ << "\" I=" << line.size() << " F=POINT" << std::endl;
                for (const auto& point : line)
                {
                    file2 << point[0] << " " << point[1] << " " << point[2] << " " << point[3] << std::endl;
                }
            }
        }
    }

    file.close();
    file2.close();
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