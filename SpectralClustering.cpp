#include "SpectralClustering.h"
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <random>
#include <limits>
#include <cmath>
#include <numeric>


SpectralClustering::SpectralClustering(const std::vector<Knot>& pts, int numClusters, float scale)
    : points(pts), k(numClusters), sigma(scale) {
}

void SpectralClustering::cluster() {
    int N = points.size();
    if (N < k) {
        throw std::invalid_argument("More clusters than points");
    }

    // 1.
    Eigen::MatrixXf S(N, N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float dx = points[i].x - points[j].x;
            float dy = points[i].y - points[j].y;
            float dist2 = dx * dx + dy * dy;
            S(i, j) = std::exp(-dist2 / (2 * sigma * sigma));
        }
    }

    // 2. 
    Eigen::VectorXf D(N);
    for (int i = 0; i < N; i++) {
        D(i) = S.row(i).sum(); 
    }

    Eigen::VectorXf D_inv_sqrt(N); 
    for (int i = 0; i < N; i++) {
        if (D(i) > 0) {
            D_inv_sqrt(i) = 1.0f / std::sqrt(D(i)); 
        }
        else {
            D_inv_sqrt(i) = 0.0f; 
        }
    }

    auto D_inv_sqrt_diag = D_inv_sqrt.asDiagonal();
    Eigen::MatrixXf S_tilde = D_inv_sqrt_diag * S * D_inv_sqrt_diag;
    Eigen::MatrixXf L = Eigen::MatrixXf::Identity(N, N) - S_tilde;

    // 3. 
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(L);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue decomposition failed");
    }

    Eigen::MatrixXf U = es.eigenvectors().block(0, 1, N, k);

    // 4. 
    for (int i = 0; i < N; i++) {
        float norm = U.row(i).norm();
        if (norm > 0) {
            U.row(i) /= norm;
        }
        else {
            U.row(i).setZero();
        }
    }

    // 5. 
    std::vector<int> indices(N);
    std::iota(indices.begin(), indices.end(), 0); 
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(indices.begin(), indices.end(), g); 

    Eigen::MatrixXf centroids(k, k); 
    for (int i = 0; i < k; i++) {
        centroids.row(i) = U.row(indices[i]);
    }

    std::vector<int> assignments(N, -1);
    int maxIter = 30;

    for (int iter = 0; iter < maxIter; iter++) {
        for (int i = 0; i < N; i++) {
            float minDist = std::numeric_limits<float>::max();
            int bestJ = -1;
            for (int j = 0; j < k; j++) {
                float dist = (U.row(i) - centroids.row(j)).squaredNorm();
                if (dist < minDist) {
                    minDist = dist;
                    bestJ = j;
                }
            }
            assignments[i] = bestJ; 
        }

        std::vector<Eigen::VectorXf> sumK(k, Eigen::VectorXf::Zero(k));
        std::vector<int> countK(k, 0); 
        for (int i = 0; i < N; i++) {
            int j = assignments[i];
            sumK[j] += U.row(i).transpose();
            countK[j]++;
        }
        for (int j = 0; j < k; j++) {
            if (countK[j] > 0) {
                sumK[j] /= countK[j]; 
                centroids.row(j) = sumK[j].transpose(); 
            }
        }
    }

    clusters = assignments; 
}

const std::vector<int>& SpectralClustering::getClusters() const {
    return clusters;
}