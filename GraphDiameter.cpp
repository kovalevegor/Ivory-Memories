#include "GraphDiameter.h"
#include <stdexcept> // Для обработки ошибок

namespace GraphAnalysis {
    // Реализация построения списка смежности
    std::unordered_map<Knot, std::vector<Knot>> buildAdjacencyList(const std::vector<Edge>& edges) {
        std::unordered_map<Knot, std::vector<Knot>> adjList;
        for (const auto& edge : edges) {
            adjList[edge.first].push_back(edge.second);
            adjList[edge.second].push_back(edge.first); // Граф неориентированный
        }
        return adjList;
    }

    // Реализация BFS
    std::unordered_map<Knot, int> bfsDistances(const Knot& start,
        const std::unordered_map<Knot, std::vector<Knot>>& adjList) {
        std::unordered_map<Knot, int> distances;
        std::queue<Knot> q;

        if (adjList.empty())
            throw std::invalid_argument("Список смежности пуст!");

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

    // Реализация поиска самой удалённой вершины
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

    // Реализация поиска диаметра графа
    std::pair<Knot, Knot> findGraphDiameter(const std::unordered_map<Knot, std::vector<Knot>>& adjList) {
        if (adjList.empty())
            return {}; // Возвращаем пустую пару для пустого графа

        // Шаг 1: Выбираем произвольную вершину
        Knot A = adjList.begin()->first;

        // Шаг 2: Находим самую удалённую вершину B от A
        auto [B, _] = findFurthestNode(A, adjList);

        // Шаг 3: Находим самую удалённую вершину C от B
        auto [C, maxDist] = findFurthestNode(B, adjList);

        return { B, C };
    }

    // Расчет эксцентриситетов для всех вершин
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

    // Поиск центра графа
    Knot findGraphCenter(const std::unordered_map<Knot, std::vector<Knot>>& adjList) {
        if (adjList.empty()) {
            throw std::invalid_argument("Граф пуст!");
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