#pragma once
#include <vector>
#include "Knot.h"

class SpectralClustering {
private:
    std::vector<Knot> points;
    int k;
    float sigma;
    std::vector<int> clusters;

public:
    SpectralClustering(const std::vector<Knot>& pts, int numClusters, float scale);
    void cluster();
    const std::vector<int>& getClusters() const;
};