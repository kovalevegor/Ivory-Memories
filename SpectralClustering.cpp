#include "SpectralClustering.h"
#include <Eigen/Dense>         // Библиотека для работы с матрицами
#include <vector>
#include <algorithm>           // Для std::shuffle
#include <random>              // Генерация случайных чисел
#include <limits>              // Для std::numeric_limits<float>::max()
#include <cmath>               // Математические функции (например, sqrt, exp)
#include <numeric>             // Для std::iota

// Конструктор: инициализирует точки, число кластеров и параметр sigma
SpectralClustering::SpectralClustering(const std::vector<Knot>& pts, int numClusters, float scale)
    : points(pts), k(numClusters), sigma(scale) {}

// Основной метод кластеризации
void SpectralClustering::cluster() {
    int N = points.size();
    if (N < k) {
        throw std::invalid_argument("More clusters than points"); // Проверка на корректность
    }

    // 1. Построение матрицы сходства S
    Eigen::MatrixXf S(N, N); // Матрица размера NxN
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float dx = points[i].x - points[j].x; // Разница по x
            float dy = points[i].y - points[j].y; // Разница по y
            float dist2 = dx * dx + dy * dy;      // Квадрат расстояния
            // Гауссово ядро: преобразует расстояние в меру сходства (0..1)
            S(i, j) = std::exp(-dist2 / (2 * sigma * sigma));
        }
    }

    // 2. Построение матрицы Лапласа L
    Eigen::VectorXf D(N); // Вектор степеней (суммы строк матрицы S)
    for (int i = 0; i < N; i++) {
        D(i) = S.row(i).sum(); // Сумма элементов в строке i
    }

    Eigen::VectorXf D_inv_sqrt(N); // Вектор обратных квадратных корней степеней
    for (int i = 0; i < N; i++) {
        if (D(i) > 0) {
            D_inv_sqrt(i) = 1.0f / std::sqrt(D(i)); // Обратный квадратный корень
        }
        else {
            D_inv_sqrt(i) = 0.0f; // Если степень нулевая
        }
    }

    // Нормализованная матрица Лапласа: L = I - D^{-1/2} * S * D^{-1/2}
    auto D_inv_sqrt_diag = D_inv_sqrt.asDiagonal(); // Преобразование вектора в диагональную матрицу
    Eigen::MatrixXf S_tilde = D_inv_sqrt_diag * S * D_inv_sqrt_diag;
    Eigen::MatrixXf L = Eigen::MatrixXf::Identity(N, N) - S_tilde;

    // 3. Вычисление собственных векторов матрицы L
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(L);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue decomposition failed"); // Ошибка разложения
    }
    // Берем первые k собственных векторов (игнорируя первый, соответствующий нулевому собственному значению)
    Eigen::MatrixXf U = es.eigenvectors().block(0, 1, N, k);

    // 4. Нормализация строк матрицы U (приведение к единичной длине)
    for (int i = 0; i < N; i++) {
        float norm = U.row(i).norm(); // Длина строки
        if (norm > 0) {
            U.row(i) /= norm; // Деление на длину для нормализации
        }
        else {
            U.row(i).setZero(); // Если строка нулевая
        }
    }

    // 5. Применение k-means к строкам матрицы U
    // Инициализация центроидов случайными точками из U
    std::vector<int> indices(N);
    std::iota(indices.begin(), indices.end(), 0); // Заполнение 0, 1, 2, ..., N-1
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(indices.begin(), indices.end(), g); // Перемешивание индексов

    Eigen::MatrixXf centroids(k, k); // Центроиды в пространстве R^k
    for (int i = 0; i < k; i++) {
        centroids.row(i) = U.row(indices[i]); // Выбор случайных строк из U
    }

    std::vector<int> assignments(N, -1); // Метки кластеров для каждой точки
    int maxIter = 30; // Максимальное число итераций k-means

    for (int iter = 0; iter < maxIter; iter++) {
        // Назначение точек ближайшим центроидам
        for (int i = 0; i < N; i++) {
            float minDist = std::numeric_limits<float>::max();
            int bestJ = -1;
            for (int j = 0; j < k; j++) {
                // Вычисление расстояния между строкой U[i] и центроидом j
                float dist = (U.row(i) - centroids.row(j)).squaredNorm();
                if (dist < minDist) {
                    minDist = dist;
                    bestJ = j;
                }
            }
            assignments[i] = bestJ; // Назначение точки i кластеру bestJ
        }

        // Обновление центроидов
        std::vector<Eigen::VectorXf> sumK(k, Eigen::VectorXf::Zero(k)); // Суммы векторов кластеров
        std::vector<int> countK(k, 0); // Число точек в кластерах
        for (int i = 0; i < N; i++) {
            int j = assignments[i];
            sumK[j] += U.row(i).transpose(); // Добавление вектора к сумме кластера
            countK[j]++;
        }
        for (int j = 0; j < k; j++) {
            if (countK[j] > 0) {
                sumK[j] /= countK[j]; // Усреднение
                centroids.row(j) = sumK[j].transpose(); // Обновление центроида
            }
        }
    }

    clusters = assignments; // Сохранение результата
}

// Возвращает метки кластеров
const std::vector<int>& SpectralClustering::getClusters() const {
    return clusters;
}