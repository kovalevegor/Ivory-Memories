#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "Knot.h"
#include "Triangulation.h"
#include "MinimumSpanningTree.h"
#include "GraphDiameter.h"
#include "SpectralClustering.h"
#include "Quest.h"

int find_knot_id(const Knot& k, const std::vector<Knot>& knots) {
    for (size_t i = 0; i < knots.size(); ++i) {
        if (k.x == knots[i].x && k.y == knots[i].y) return i;
    }
    return -1;
}

int main() {
    srand(time(0));
    std::random_device rd;
    std::mt19937 rng(rd());

    // Параметры генерации графа
    int n = 100;                // Number of knots
    float width = 192.0f;       // Maximum x-coordinate
    float height = 108.0f;      // Maximum y-coordinate
    float radius = 10;          // Minimum radius between knots
    int k = 30;             // Maximum attempts 

    std::vector<Knot> knots;
    knots.reserve(n);

    knots = poissonDiskSampling(width, height, radius, k);
    for (const auto& knot : knots) {
        std::cout << "(" << knot.x << ", " << knot.y << ")\n";
    }

    {
        std::ofstream out("knots.txt");
        for (const auto& knot : knots) {
            out << knot.x << " " << knot.y << "\n";
        }
    }

    std::vector<Triangle> triangles = delaunayTriangulation(knots);
    std::cout << "Number of triangles: " << triangles.size() << std::endl;
    for (const auto& triangle : triangles) {
        std::cout << "Triangle: ("
            << triangle.knot_1.x << ", " << triangle.knot_1.y << "), ("
            << triangle.knot_2.x << ", " << triangle.knot_2.y << "), ("
            << triangle.knot_3.x << ", " << triangle.knot_3.y << ")\n";
    }

    {
        std::ofstream out("triangles.txt");
        for (const auto& triangle : triangles) {
            out << triangle.knot_1.x << " " << triangle.knot_1.y << " "
                << triangle.knot_2.x << " " << triangle.knot_2.y << " "
                << triangle.knot_3.x << " " << triangle.knot_3.y << "\n";
        }
    }

    std::vector<Edge> delaunay_polygon = getUniqueEdges(triangles);
    {
        std::ofstream out("delaunay_polygon.txt");
        for (const auto& edge : delaunay_polygon) {
            out << edge.first.x << " " << edge.first.y << " "
                << edge.second.x << " " << edge.second.y << "\n";
        }
    }

    std::vector<EdgeWithWeight> mst = calculateMST(delaunay_polygon);
    {
        std::ofstream out("mst.txt");
        for (const auto& edge : mst) {
            out << edge.start.x << " " << edge.start.y << " "
                << edge.end.x << " " << edge.end.y << "\n";
        }
    }
    // Преобразование EdgeWithWeight в Edge
    std::vector<Edge> mst_edges;
    for (const auto& weighted_edge : mst) {
        mst_edges.emplace_back(weighted_edge.start, weighted_edge.end);
    }

    float probability = 0.25f;
    auto maze = addRandomEdges(mst_edges, delaunay_polygon, probability);
    {
        std::ofstream out("maze.txt");
        for (const auto& edge : maze) {
            out << edge.first.x << " " << edge.first.y << " "
                << edge.second.x << " " << edge.second.y << "\n";
        }
    }

    // Построение списка смежности
    auto adjacencyList = GraphAnalysis::buildAdjacencyList(maze);

    // Поиск диаметра графа
    auto [spawn, exit] = GraphAnalysis::findGraphDiameter(adjacencyList);

    std::cout << "Spawn knot: (" << spawn.x << ", " << spawn.y << ")\n";
    std::cout << "Exit knot: (" << exit.x << ", " << exit.y << ")\n";
    {
        std::ofstream out("spawn_exit_knots.txt");
        out << spawn.x << " " << spawn.y << "\n"
            << exit.x << " " << exit.y << "\n";
    }

    // Поиск центра графа (убежища)
    try {
        Knot center = GraphAnalysis::findGraphCenter(adjacencyList);
        std::cout << "Graph center: (" << center.x << ", " << center.y << ")\n";
        {
            std::ofstream out("asylum.txt");
            out << center.x << " " << center.y << "\n";
        }
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        {
            std::ofstream out("asylum.txt");
            out << "Error: " << e.what() << "\n";
        }
    }

    // Кластеризация
    SpectralClustering clustering(knots, 5, 10.0f);
    std::vector<int> cluster_assignments;

    try {
        clustering.cluster();
        cluster_assignments = clustering.getClusters();

        std::ofstream clustering_file("clustering.txt");
        for (size_t i = 0; i < knots.size(); ++i) {
            std::string line = "Knot (" + std::to_string(knots[i].x) + ", " + std::to_string(knots[i].y) +
                ") -> Cluster " + std::to_string(cluster_assignments[i]);
            std::cout << line << std::endl;
            clustering_file << line << "\n";
        }
        clustering_file.close();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::ofstream clustering_file("clustering.txt");
        clustering_file << "Error: " << e.what() << "\n";
        clustering_file.close();
    }


    // Распределение квестов
    // Создание игрового мира
    World world;
    world.knots = knots;

    // Конвертация списка смежности в формат с ID узлов
    for (const auto& [knot, neighbors] : adjacencyList) {
        int knot_id = find_knot_id(knot, knots);
        std::vector<int> neighbor_ids;
        for (const auto& neighbor : neighbors) {
            int neighbor_id = find_knot_id(neighbor, knots);
            neighbor_ids.push_back(neighbor_id);
        }
        world.adjacency_list[knot_id] = neighbor_ids;
    }

    // Размещение объектов по кластерам
    const int num_items = 5;
    const int num_characters = 5; // По одному персонажу на кластер
    const int num_enemies = 5;    // По одному врагу на кластер

    // Распределение узлов по кластерам
    std::vector<std::vector<int>> clusters(5);
    for (size_t i = 0; i < knots.size(); ++i) {
        if (cluster_assignments[i] >= 0 && cluster_assignments[i] < 5) {
            clusters[cluster_assignments[i]].push_back(i);
        }
    }

    // Размещение персонажей (по одному в каждом кластере)
    for (int i = 0; i < num_characters; ++i) {
        if (!clusters[i].empty()) {
            int knot_id = clusters[i][rand() % clusters[i].size()];
            world.character_locations[i] = knot_id;
            // Случайное требование предмета или отсутствие требования
            world.character_requirements[i] = (rand() % 2 == 0) ? (rand() % num_items) : -1;
        }
    }

    // Размещение врагов (по одному в каждом кластере)
    for (int i = 0; i < num_enemies; ++i) {
        if (!clusters[i].empty()) {
            int knot_id = clusters[i][rand() % clusters[i].size()];
            world.enemy_locations[i] = knot_id;
        }
    }

    // Размещение источников предметов
    for (int item = 0; item < num_items; ++item) {
        std::vector<Source> sources;
        // Источник на узле (в случайном кластере)
        int cluster_id = rand() % 5;
        if (!clusters[cluster_id].empty()) {
            int knot_id = clusters[cluster_id][rand() % clusters[cluster_id].size()];
            sources.push_back({ NODE_SOURCE, knot_id });
        }
        // Источник от персонажа (случайный персонаж)
        int char_id = rand() % num_characters;
        sources.push_back({ CHARACTER_SOURCE, char_id });
        // Источник от врага (случайный враг)
        int enemy_id = rand() % num_enemies;
        sources.push_back({ ENEMY_SOURCE, enemy_id });
        world.item_sources[item] = sources;
    }

    // Определение квеста
    int quest_giver_knot = find_knot_id(spawn, knots);
    int fetch_item = 0; // Цель — добыть предмет 0
    State initial;
    initial.player_location = quest_giver_knot;
    for (int item = 0; item < num_items; ++item) {
        initial.item_locations[item] = world.item_sources[item][0].id; // Начальное расположение в первом источнике
    }

    // Запуск генетического алгоритма
    Individual best_individual = genetic_algorithm(50, 100, world, fetch_item, quest_giver_knot, initial, rng);

    // Построение и сохранение дерева квестов
    std::unordered_map<std::pair<State, Goal>, std::vector<Action>> plan_cache;
    QuestNode* quest_tree = build_fetch_quest_tree(fetch_item, quest_giver_knot, best_individual, world, initial, plan_cache);
    save_quest_tree(quest_tree, "quest_tree.txt");
    delete quest_tree;

    return 0;
}