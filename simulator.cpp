// ==========================
// TITLE: DCFHT PIC SIMULATOR (Optimized)
// AUTHOR: ZACHARY LINDER
// VERSION 4.5
// RECENT DATE: 3.9.2025
// ==========================
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <cstdlib>
#include <sstream>
#include "progressbar.hpp"
#include <omp.h>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <execution>
#include <chrono>
#include <iomanip>

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

// -----------------------------------------------------------------------------
// Precomputed constants to avoid repeated divisions
// -----------------------------------------------------------------------------
constexpr double CELL_RESOLUTION = 0.0005;
constexpr double INV_CELL_RESOLUTION = 1.0 / CELL_RESOLUTION;  // 1000.0
constexpr double DIVIDE_BY_10 = 0.1;

// ===============================
// === RANDOM NUMBER GENERATOR ===
// ===============================
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis;
inline double random_double(double min, double max)
{
    return dis(gen, std::uniform_real_distribution<double>::param_type(min, max));
}

// ==========================
// === STRUCTURED VECTORS ===
// ==========================
struct Cell
{
    int A, B, C, D; // Represent 4 Indexes Per Cell
    double Aw, Bw, Cw, Dw; // Represent 4 Weights Per Cell
};

struct Vector3
{
    double x, y, z; // 3-element vector used for positions and velocities
};

struct Field
{
    double x, y, z, mx, my, mz; // Field data parameters
};

struct ReflectionData
{
    double normalX, normalY; // Wall-normal vector components
};

// =======================================================
// === DETERMINE SIZE OF FIELD DATA (ROWS AND COLUMNS) ===
// =======================================================
bool compareFields(const Field& a, const Field& b) {
    return (a.x != b.x) ? a.x < b.x : a.y < b.y;
}

int countYPerX(std::vector<Field>& bField) {
    std::sort(bField.begin(), bField.end(), compareFields);
    int yCount = 0;
    // Count unique y values for the first x value only.
    for (size_t i = 0; i < bField.size();) {
        double currentX = bField[i].x;
        double lastY = bField[i].y - 1; // ensure difference
        while (i < bField.size() && bField[i].x == currentX) {
            if (bField[i].y != lastY) {
                ++yCount;
                lastY = bField[i].y;
            }
            ++i;
        }
        break; // only count for the first x value
    }
    return yCount;
}

// ============================
// === READ FIELD FILE DATA ===
// ============================
std::vector<Field> read_fld_file(const std::string& filename)
{
    std::ifstream file(filename);
    std::vector<Field> field;
    field.reserve(10000);
    std::string line;
    int lineCount = 0;
    while (std::getline(file, line))
    {
        if (++lineCount < 2) continue; // Skip header lines
        Field point;
        std::istringstream iss(line);
        if (iss >> point.x >> point.y >> point.z >> point.mx >> point.my >> point.mz) {
            field.push_back(point);
        }
    }
    return field;
}

// ==============================
// === MATRIX MATH OPERATIONS ===
// ==============================
void MCM(std::vector<Vector3>& result, const std::vector<Vector3>& matrix, double constant)
{
    std::transform(std::execution::par, matrix.begin(), matrix.end(), result.begin(),
        [constant](const Vector3& v) { return Vector3{ v.x * constant, v.y * constant, v.z * constant }; });
}

void MCMA(std::vector<Vector3>& result, const std::vector<Vector3>& matrix, const std::vector<double>& constant)
{
    std::transform(std::execution::par, matrix.begin(), matrix.end(), constant.begin(), result.begin(),
        [](const Vector3& vec, double scalar) { return Vector3{ vec.x * scalar, vec.y * scalar, vec.z * scalar }; });
}

void MCA(std::vector<Vector3>& result, const std::vector<Vector3>& matrix, const std::vector<double>& constant)
{
    std::transform(std::execution::par, matrix.begin(), matrix.end(), constant.begin(), result.begin(),
        [](const Vector3& vec, double scalar) { return Vector3{ vec.x + scalar, vec.y + scalar, vec.z + scalar }; });
}

void CCAA(std::vector<double>& result, const std::vector<double>& C1, double C2)
{
    std::transform(std::execution::par, C1.begin(), C1.begin() + 100, result.begin(),
        [C2](double val) { return val + C2; });
}

void CCDA(std::vector<double>& result, double C1, const std::vector<double>& C2)
{
    std::transform(std::execution::par, C2.begin(), C2.begin() + 100, result.begin(),
        [C1](double val) { return C1 / val; });
}

void MMA(std::vector<Vector3>& result, const std::vector<Vector3>& matrixA, const std::vector<Vector3>& matrixB)
{
    std::transform(std::execution::par, matrixA.begin(), matrixA.end(), matrixB.begin(), result.begin(),
        [](const Vector3& a, const Vector3& b) { return Vector3{ a.x + b.x, a.y + b.y, a.z + b.z }; });
}

void MMS(std::vector<double>& result, const std::vector<Vector3>& matrix)
{
    result.resize(matrix.size());
    std::transform(std::execution::par, matrix.begin(), matrix.end(), result.begin(),
        [](const Vector3& vec) { return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z; });
}

void CROSS(std::vector<Vector3>& result, const std::vector<Vector3>& matrixA, const std::vector<Vector3>& matrixB)
{
    result.resize(matrixA.size());
    std::transform(std::execution::par, matrixA.begin(), matrixA.end(), matrixB.begin(), result.begin(),
        [](const Vector3& a, const Vector3& b) -> Vector3 {
            return { a.y * b.z - a.z * b.y,
                     -(a.x * b.z - a.z * b.x),
                     a.x * b.y - a.y * b.x };
        });
}

// ======================
// === WALL COLLISION ===
// ======================
inline bool isOutOfBounds(const Vector3& position, double slope_left, double slope_right,
    double y_intercept_left, double y_intercept_right,
    double Y1, double Y2, int& type, double center)
{
    if (position.y > Y1) {
        type = 4;
        return true;
    }
    else if (position.y < Y2) {
        type = 2;
        return true;
    }
    else if (position.x < (position.y - y_intercept_left) / slope_left && position.x < center) {
        type = 1;
        return true;
    }
    else if (position.x > (position.y - y_intercept_right) / slope_right && position.x > center) {
        type = 3;
        return true;
    }
    return false;
}

// ===================================
// === REFLECTION VIA WALL NORMALS ===
// ===================================
inline Vector3 reflect(const Vector3& vec, double normalX, double normalY)
{
    double invNorm = 1.0 / std::sqrt(normalX * normalX + normalY * normalY);
    normalX *= invNorm;
    normalY *= invNorm;
    double dotProduct = vec.x * normalX + vec.y * normalY;
    double scale = 2 * dotProduct;
    return { vec.x - scale * normalX, vec.y - scale * normalY, vec.z };
}

// =========================
// === POSITION --> CELL ===
// =========================
inline int findCell(const Vector3& particle, double dx, double dy, double minX, double minY, int numX, int numY)
{
    // Instead of using floor we use multiplication by precomputed inverse.
    int col = static_cast<int>((particle.x - minX) * INV_CELL_RESOLUTION);
    int row = static_cast<int>((particle.y - minY) * INV_CELL_RESOLUTION);
    const int numRows = numX - 1;
    return (col - 1) * numRows + row;
}

// ==========================
// === FIELD --> PARTICLE ===
// ==========================
void FieldSolve(std::vector<Vector3>& particlePos, std::vector<Vector3>& gridPos,
    std::vector<Vector3>& eMag, std::vector<Vector3>& bMag,
    std::vector<Vector3>& E, std::vector<Vector3>& B,
    int cells[], std::vector<Cell>& grid,
    double dx, double dy, double minX, double minY,
    int discrete, int t, int numX, int numY)
{
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(particlePos.size()); ++i)
    {
        int cell = findCell(particlePos[i], dx, dy, minX, minY, numX, numY);
        // Use cell-1 for 0-based index update.
#pragma omp atomic
        cells[cell - 1] += 1;

        // Precompute differences for the four surrounding grid points.
        int Ag = grid[cell].A, Bg = grid[cell].B, Cg = grid[cell].C, Dg = grid[cell].D;
        auto calcRatio = [&](int idx) -> long double {
            double diffX = gridPos[idx].x - particlePos[i].x;
            double diffY = gridPos[idx].y - particlePos[i].y;
            // Avoid division by zero.
            return (std::abs(diffX) < 1e-12) ? 0.0L : (diffY / diffX);
            };
        long double A1 = calcRatio(Ag);
        long double A2 = calcRatio(Bg);
        long double A3 = calcRatio(Cg);
        long double A4 = calcRatio(Dg);

        long double mag = std::sqrt(A1 * A1 + A2 * A2 + A3 * A3 + A4 * A4);
        if (mag < 1e-10L) mag = 1e-10L;
        Cell weight = { 0, 0, 0, 0,
            static_cast<double>(std::abs(A1 / mag)),
            static_cast<double>(std::abs(A2 / mag)),
            static_cast<double>(std::abs(A3 / mag)),
            static_cast<double>(std::abs(A4 / mag)) };

        //E[i] = { eMag[Ag].x * weight.Aw + eMag[Bg].x * weight.Bw +
                 //eMag[Cg].x * weight.Cw + eMag[Dg].x * weight.Dw,
                 //eMag[Ag].y * weight.Aw + eMag[Bg].y * weight.Bw +
                 //eMag[Cg].y * weight.Cw + eMag[Dg].y * weight.Dw,
                 //0 };
        E[i] = { -10,-800,0 };
        B[i] = { bMag[Ag].x * weight.Aw + bMag[Bg].x * weight.Bw +
                 bMag[Cg].x * weight.Cw + bMag[Dg].x * weight.Dw,
                 bMag[Ag].y * weight.Aw + bMag[Bg].y * weight.Bw +
                 bMag[Cg].y * weight.Cw + bMag[Dg].y * weight.Dw,
                 0 };
    }
    if (t % discrete == 0)
    {
        std::ofstream outFile("discrete.txt", std::ios::app);
        if (outFile.is_open())
        {
            for (int j = 0; j < numY - 1; j++)
            {
                for (int i = 0; i < numX - 1; ++i)
                    outFile << cells[j * (numX - 1) + i] << " ";
                outFile << "\n";
            }
            outFile << "\n";
            outFile.close();
        }
    }
}

// ==========================
// === PARTICLE --> FIELD ===
// ==========================
void ParticleSolve(std::vector<Vector3>& particlePos, std::vector<Vector3>& gridPos,
    std::vector<Vector3>& eMag, std::vector<Vector3>& bMag,
    std::vector<Vector3>& E, std::vector<Vector3>& B,
    int cells[], std::vector<Cell>& grid,
    double dx, double dy, double minX, double minY,
    int discrete, int t, int numX, int numY)
{
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(particlePos.size()); ++i)
    {
        int cell = findCell(particlePos[i], dx, dy, minX, minY, numX, numY);
        int Ag = grid[cell].A, Bg = grid[cell].B, Cg = grid[cell].C, Dg = grid[cell].D;
        auto calcRatio = [&](int idx) -> long double {
            double diffX = gridPos[idx].x - particlePos[i].x;
            double diffY = gridPos[idx].y - particlePos[i].y;
            return (std::abs(diffX) < 1e-12) ? 0.0L : (diffY / diffX);
            };
        long double A1 = calcRatio(Ag);
        long double A2 = calcRatio(Bg);
        long double A3 = calcRatio(Cg);
        long double A4 = calcRatio(Dg);

        long double mag = std::sqrt(A1 * A1 + A2 * A2 + A3 * A3 + A4 * A4);
        if (mag < 1e-10L) mag = 1e-10L;
        Cell weight = { 0, 0, 0, 0,
            static_cast<double>(std::abs(A1 / mag)),
            static_cast<double>(std::abs(A2 / mag)),
            static_cast<double>(std::abs(A3 / mag)),
            static_cast<double>(std::abs(A4 / mag)) };

        eMag[Ag].x += weight.Aw * E[i].x * DIVIDE_BY_10;
        eMag[Ag].y += weight.Aw * E[i].y * DIVIDE_BY_10;
        eMag[Bg].x += weight.Bw * E[i].x * DIVIDE_BY_10;
        eMag[Bg].y += weight.Bw * E[i].y * DIVIDE_BY_10;
        eMag[Cg].x += weight.Cw * E[i].x * DIVIDE_BY_10;
        eMag[Cg].y += weight.Cw * E[i].y * DIVIDE_BY_10;
        eMag[Dg].x += weight.Dw * E[i].x * DIVIDE_BY_10;
        eMag[Dg].y += weight.Dw * E[i].y * DIVIDE_BY_10;
        bMag[Ag].x += weight.Aw * B[i].x * DIVIDE_BY_10;
        bMag[Ag].y += weight.Aw * B[i].y * DIVIDE_BY_10;
        bMag[Bg].x += weight.Bw * B[i].x * DIVIDE_BY_10;
        bMag[Bg].y += weight.Bw * B[i].y * DIVIDE_BY_10;
        bMag[Cg].x += weight.Cw * B[i].x * DIVIDE_BY_10;
        bMag[Cg].y += weight.Cw * B[i].y * DIVIDE_BY_10;
        bMag[Dg].x += weight.Dw * B[i].x * DIVIDE_BY_10;
        bMag[Dg].y += weight.Dw * B[i].y * DIVIDE_BY_10;
    }
}

// =============================
// === GRID POINTS --> CELLS ===
// =============================
std::vector<Cell> generateGrid(int numX, int numY)
{
    int num_y = numY; // number of x values per row
    int num_x = numX; // number of y values per column
    std::vector<Cell> grid((num_x - 1) * (num_y - 1) + 1);
    for (int i = 0; i < num_y - 1; ++i) {
        for (int j = 0; j < num_x - 1; ++j) {
            int base = j + i * num_x;
            Cell cell;
            cell.A = base + 1;
            cell.B = base + 2;
            cell.C = base + num_x + 1;
            cell.D = base + num_x + 2;
            grid[j + i * (num_x - 1) + 1] = cell; // 1-based indexing
        }
    }
    return grid;
}

// =====================
// === MAIN FUNCTION ===
// =====================
int main()
{
    auto start_time = std::chrono::high_resolution_clock::now();

    // === Read field data ===
    std::vector<Field> bField = read_fld_file("bField.fld");
    std::vector<Field> eField = read_fld_file("eField.fld");

    // === Read settings ===
    std::ifstream settingsFile("settings.txt");
    std::unordered_map<std::string, double> settings;
    std::string line;
    while (std::getline(settingsFile, line))
    {
        std::istringstream iss(line);
        std::string key;
        double value;
        if (std::getline(iss, key, '=') && (iss >> value))
        {
            key.erase(key.find_last_not_of(" \t") + 1);
            settings[key] = value;
        }
    }
    settingsFile.close();

    int particles = static_cast<int>(settings["PARTICLES"]);
    int timesteps = static_cast<int>(settings["TIMESTEPS"]);
    int points = static_cast<int>(settings["POINTS"]);
    double dt = settings["DT"];
    double qom = settings["QOM"];
    double outletRadius = settings["RADIUS"];
    double div = settings["DIV"];
    double inletRadius = settings["INLET"];
    int plot = settings["PLOT"];

    std::cout << "DCFHT PIC SIM" << std::endl;
    std::cout << "Particles: " << particles << std::endl;
    std::cout << "Simulation Time (Absolute): " << timesteps * dt << " s" << std::endl;

    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    if (!bField.empty()) {
        auto [minXIter, maxXIter] = std::minmax_element(bField.begin(), bField.end(),
            [](const Field& a, const Field& b) { return a.x < b.x; });
        auto [minYIter, maxYIter] = std::minmax_element(bField.begin(), bField.end(),
            [](const Field& a, const Field& b) { return a.y < b.y; });
        minX = minXIter->x;
        maxX = maxXIter->x;
        minY = minYIter->y;
        maxY = maxYIter->y;
    }

    double divRad = div * M_PI / 180;
    double center = (minX + maxX) * 0.5;
    double outletMinX = center - outletRadius;
    double outletMaxX = center + outletRadius;
    double outletY = maxY;
    double inletMinX = center - inletRadius;
    double inletMaxX = center + inletRadius;
    double NormalLeftX = std::cos(divRad * 0.5);
    double NormalLeftY = std::sin(divRad * 0.5);
    double NormalRightX = -NormalLeftX;
    double NormalRightY = NormalLeftY;
    double inletY = outletY - ((inletMinX - outletMinX) * std::sin(M_PI * 0.5 - divRad * 0.5) / std::sin(divRad * 0.5));

    std::cout << "THRUSTER GEOMETRY:" << std::endl;
    std::cout << "\tHeight: " << (outletY - inletY) * 1e3 << " mm" << std::endl;
    std::cout << "\tOutlet Width: " << 2 * outletRadius * 1e3 << " mm" << std::endl;
    std::cout << "\tInlet Width: " << 2 * inletRadius * 1e3 << " mm" << std::endl;
    std::cout << "\tDivergence Angle: " << div << " deg" << std::endl;

    double slope_left = (outletY - inletY) / (outletMinX - inletMinX);
    double slope_right = (outletY - inletY) / (outletMaxX - inletMaxX);
    double y_intercept_left = outletY - (slope_left * outletMinX);
    double y_intercept_right = outletY - (slope_right * outletMaxX);

    // === Allocate simulation arrays ===
    std::vector<Vector3> particlePos(particles, { 0.0, 0.0, 0.0 });
    std::vector<Vector3> particleVel(particles, { 0.0, 0.0, 0.0 });
    std::vector<Vector3> bPos(bField.size()), ePos(eField.size()),
        bMag(bField.size()), eMag(bField.size());
    std::vector<Vector3> E(particles, { 0.0, -300, 0.0 }), B(particles, { 0.0, -1, 0.0 });
    std::vector<Vector3> vm(particles), t(particles), vd(particles), vp(particles),
        E_term(particles), vel_term(particles), s(particles),
        vd_term(particles), vm_term(particles);
    std::vector<double> t_term3(particles, 0.0), t_term2(particles, 0.0), t_term(particles, 0.0);

    std::unordered_map<int, ReflectionData> reflectionDataMap =
    {
        {0, {0.0, 0.0}},
        {1, {NormalLeftX, NormalLeftY}},
        {2, {0, 1}},
        {3, {NormalRightX, NormalRightY}}
    };

    // === Initialize particle positions and velocities in parallel ===
    #pragma omp parallel for
    for (int i = 0; i < particles; ++i)
    {
        particlePos[i] = { random_double(outletMinX + 0.001, outletMaxX - 0.001),
                           outletY - 0.001, 0.0 };
        particleVel[i] = { 0.0, 0.0, 0.0 };
    }
    #pragma omp parallel for
    for (size_t i = 0; i < bField.size(); ++i)
    {
        ePos[i] = { eField[i].x, eField[i].y, 0 };
        bPos[i] = { bField[i].x, bField[i].y, 0 };
        eMag[i] = { eField[i].mx, eField[i].my, 0 };
        bMag[i] = { bField[i].mx, bField[i].my, 0 };
    }

    std::ofstream logFile("log.txt", std::ios::trunc);
    std::ofstream outFile("discrete.txt", std::ios::trunc);
    int numX = countYPerX(bField);
    int numY = static_cast<int>(bField.size()) / numX;
    std::cout << "Grid Size: " << numX << " X " << numY << std::endl;
    double dx = (minX - bPos[numX * numY].x) / numX;
    double dy = (minY - bPos[numX * numY].y) / numY;
    int cells[(numX - 1) * (numY - 1)] = {};
    int recycled = 0;
    int precision = timesteps / points;
    double t_t = qom * dt * 0.5;
    int discrete = timesteps / plot;

    std::vector<Cell> grid = generateGrid(numX, numY);
    progressbar bar(100);

    // === Main simulation loop ===
    for (int i = 0; i < timesteps; ++i)
    {
#pragma omp parallel for
        for (size_t a = 0; a < bField.size(); ++a)
        {
            eMag[a] = { eField[a].mx, eField[a].my, 0 };
            bMag[a] = { bField[a].mx, bField[a].my, 0 };
        }

        ParticleSolve(particlePos, bPos, eMag, bMag, E, B, cells, grid, dx, dy, minX, minY, discrete, i, numX, numY);
        FieldSolve(particlePos, bPos, eMag, bMag, E, B, cells, grid, dx, dy, minX, minY, discrete, i, numX, numY);
        MCM(E_term, E, t_t);
        MMA(vm, particleVel, E_term);
        MCM(t, B, t_t);
        MMS(t_term3, t);
        CCAA(t_term2, t_term3, 1.0);
        CCDA(t_term, 2.0, t_term2);
        MCMA(s, t, t_term);
        CROSS(vm_term, vm, t);
        MMA(vd, vm, vm_term);
        CROSS(vd_term, vd, s);
        MMA(vp, vm, vd_term);

#pragma omp parallel for
        for (int p = 0; p < particles; ++p)
        {
            int type = 0;
            if (isOutOfBounds(particlePos[p], slope_left, slope_right, y_intercept_left, y_intercept_right, maxY, inletY, type, center))
            {
                if (type == 1 || type == 2 || type == 3)
                {
                    auto& reflection = reflectionDataMap[type];
                    particleVel[p] = reflect(particleVel[p], reflection.normalX, reflection.normalY);
                }
                else if (type == 4)
                {
                    // Recycle particle
                    particlePos[p] = { random_double(outletMinX + 0.001, outletMaxX - 0.001),
                                       outletY - 0.001, 0.0 };
                    particleVel[p] = { 0.0, 0.0, 0.0 };
                    E[p] = { 0.0, 0.0, 0.0 };
                    B[p] = { 0.0, 0.0, 0.0 };
    #pragma omp atomic
                    recycled++;
                }
            }
            else {
                particleVel[p] = vp[p];
            }
        }
        MCM(vel_term, particleVel, dt);
        MMA(particlePos, particlePos, vel_term);

        if (i % precision == 0)
        {
            std::ostringstream logBuffer;
            for (int p = 0; p < particles; ++p)
            {
                logBuffer << p + 1 << " " << i / precision + 1 << " "
                    << particlePos[p].x << " " << particlePos[p].y << " " << particlePos[p].z << "\n";
            }
            logFile << logBuffer.str();
        }
        if (i % (timesteps / 100) == 0)
            bar.update();
    }

    int sum = 0;
#pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < (numX - 1) * (numY - 1); ++i)
        sum += cells[i];

    logFile.close();
    outFile.close();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    std::cout << "\nConfirmed Particle Count: " << sum / timesteps << std::endl;
    std::cout << "Number of Particles Recycled: " << recycled << std::endl;
    std::cout << "Execution time: " << duration << " ms" << std::endl;
    std::cout << "Approximate time per particle: " << duration / particles << " ms" << std::endl;
    std::cout << "Approximate time per million particles: "
        << (duration / particles) * 1e6 / 3600 / 1000 << " hrs" << std::endl;
    return 0;
}
