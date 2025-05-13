#include "GraphDiameter.h"
#include <stdexcept> 

namespace GraphAnalysis {
    std::unordered_map<Knot, std::vector<Knot>> buildAdjacencyList(const std::vector<Edge>& edges) {
        std::unordered_map<Knot, std::vector<Knot>> adjList;
        for (const auto& edge : edges) {
            adjList[edge.first].push_back(edge.second);
            adjList[edge.second].push_back(edge.first);
        }
        return adjList;
    }

    std::unordered_map<Knot, int> bfsDistances(const Knot& start,
        const std::unordered_map<Knot, std::vector<Knot>>& adjList) {
        std::unordered_map<Knot, int> distances;
        std::queue<Knot> q;

        if (adjList.empty())
            throw std::invalid_argument("Adjacency list is empty!!");

        distances[start] = 0;
        q.push(start);

        while (!q.empty()) {
            Knot current = q.front();
            q.pop();

            for (const auto& neighbor : adjList.at(current)) {
                if (!distances.contains(neighbor)) {
                    distances[neighbor] = distances[current] + 1;
                    q.push(neighbor);
                }
            }
        }

        return distances;
    }

    std::pair<Knot, int> findFurthestNode(const Knot& start,
        const std::unordered_map<Knot, std::vector<Knot>>& adjList) {
        auto distances = bfsDistances(start, adjList);
        Knot furthestNode = start;
        int maxDist = 0;

        for (const auto& [node, dist] : distances) {
            if (dist > maxDist) {
                maxDist = dist;
                furthestNode = node;
            }
        }

        return { furthestNode, maxDist };
    }

    std::pair<Knot, Knot> findGraphDiameter(const std::unordered_map<Knot, std::vector<Knot>>& adjList) {
        if (adjList.empty())
            return {}; 

        Knot A = adjList.begin()->first;
        auto [B, _] = findFurthestNode(A, adjList);
        auto [C, maxDist] = findFurthestNode(B, adjList);

        return { B, C };
    }

    std::unordered_map<Knot, int> calculateEccentricities(const std::unordered_map<Knot, std::vector<Knot>>& adjList) {
        std::unordered_map<Knot, int> eccentricities;

        for (const auto& [node, _] : adjList) {
            auto distances = bfsDistances(node, adjList);
            int maxDist = 0;
            for (const auto& [_, dist] : distances) {
                if (dist > maxDist) maxDist = dist;
            }
            eccentricities[node] = maxDist;
        }

        return eccentricities;
    }

    Knot findGraphCenter(const std::unordered_map<Knot, std::vector<Knot>>& adjList) {
        if (adjList.empty()) {
            throw std::invalid_argument("Adjacency list is empty!");
        }

        auto eccentricities = calculateEccentricities(adjList);
        Knot center;
        int minEccentricity = std::numeric_limits<int>::max();

        for (const auto& [node, ecc] : eccentricities) {
            if (ecc < minEccentricity) {
                minEccentricity = ecc;
                center = node;
            }
        }

        return center;
    }
}