#pragma once

#pragma once
#include "Knot.h"
#include "Triangulation.h"
#include <vector>
#include <unordered_map>
#include <queue>
#include <utility> // Для std::pair

// Структура для хеширования Edge (если требуется)
struct EdgeHash;

// Основные функции для работы с графом
namespace GraphAnalysis {
    // Преобразование списка рёбер в список смежности
    std::unordered_map<Knot, std::vector<Knot>> buildAdjacencyList(const std::vector<Edge>& edges);

    // Обход в ширину (BFS) для вычисления расстояний
    std::unordered_map<Knot, int> bfsDistances(const Knot& start,
        const std::unordered_map<Knot, std::vector<Knot>>& adjList);

    // Поиск самой удалённой вершины от заданной
    std::pair<Knot, int> findFurthestNode(const Knot& start,
        const std::unordered_map<Knot, std::vector<Knot>>& adjList);

    // Поиск диаметра графа (пары самых удалённых вершин)
    std::pair<Knot, Knot> findGraphDiameter(const std::unordered_map<Knot, std::vector<Knot>>& adjList);

    // Найти эксцентриситет для каждой вершины (максимальное расстояние до других вершин)
    std::unordered_map<Knot, int> calculateEccentricities(const std::unordered_map<Knot, std::vector<Knot>>& adjList);

    // Найти центр графа (вершину с минимальным эксцентриситетом)
    Knot findGraphCenter(const std::unordered_map<Knot, std::vector<Knot>>& adjList);
}