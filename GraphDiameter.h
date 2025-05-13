#pragma once

#pragma once
#include "Knot.h"
#include "Triangulation.h"
#include <vector>
#include <unordered_map>
#include <queue>
#include <utility> 

struct EdgeHash;

namespace GraphAnalysis {
    std::unordered_map<Knot, std::vector<Knot>> buildAdjacencyList(const std::vector<Edge>& edges);

    std::unordered_map<Knot, int> bfsDistances(const Knot& start,
        const std::unordered_map<Knot, std::vector<Knot>>& adjList);

    std::pair<Knot, int> findFurthestNode(const Knot& start,
        const std::unordered_map<Knot, std::vector<Knot>>& adjList);

    std::pair<Knot, Knot> findGraphDiameter(const std::unordered_map<Knot, std::vector<Knot>>& adjList);

    std::unordered_map<Knot, int> calculateEccentricities(const std::unordered_map<Knot, std::vector<Knot>>& adjList);

    Knot findGraphCenter(const std::unordered_map<Knot, std::vector<Knot>>& adjList);
}