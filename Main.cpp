#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "Knot.h"
#include "Triangulation.h"
#include "MinimumSpanningTree.h"
#include "GraphDiameter.h"
#include "SpectralClustering.h"
#include "QuestPlacement.h"

#include <json.hpp>

using json = nlohmann::json;

int find_knot_id(const Knot& k, const std::vector<Knot>& knots) {
    const float EPSILON = 1e-5;
    for (size_t i = 0; i < knots.size(); ++i) {
        if (std::abs(k.x - knots[i].x) < EPSILON && std::abs(k.y - knots[i].y) < EPSILON) return i;
    }
    return -1;
}

int main() {
    srand(time(0));
    std::random_device rd;
    std::mt19937 rng(rd());

    // Параметры генерации графа
    //int n = 100;                // Number of knots
    float width = 192.0f;       // Maximum x-coordinate
    float height = 108.0f;      // Maximum y-coordinate
    float radius = 10;          // Minimum radius between knots
    int k = 30;                 // Maximum attempts 

    std::vector<Knot> knots;
    //knots.reserve(n);

    knots = poissonDiskSampling(width, height, radius, k);
    for (const auto& knot : knots) {
        //std::cout << "(" << knot.x << ", " << knot.y << ")\n";
    }
    std::cout << "Number of knots: " << knots.size() << std::endl;

    {
        std::ofstream out("knots.txt");
        for (const auto& knot : knots) {
            out << knot.x << " " << knot.y << "\n";
        }
    }

    std::vector<Triangle> triangles = delaunayTriangulation(knots);
    std::cout << "Number of triangles: " << triangles.size() << std::endl;
    /*for (const auto& triangle : triangles) {
        std::cout << "Triangle: ("
            << triangle.knot_1.x << ", " << triangle.knot_1.y << "), ("
            << triangle.knot_2.x << ", " << triangle.knot_2.y << "), ("
            << triangle.knot_3.x << ", " << triangle.knot_3.y << ")\n";
    }*/

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
    Knot center;
    try {
        center = GraphAnalysis::findGraphCenter(adjacencyList);
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

    // Проверка, что spawn knot присутствует в maze 
    bool spawn_in_maze = false;
    for (const auto& edge : maze) {
        if ((edge.first.x == spawn.x && edge.first.y == spawn.y) ||
            (edge.second.x == spawn.x && edge.second.y == spawn.y)) {
            spawn_in_maze = true;
            break;
        }
    }
    if (!spawn_in_maze) {
        std::cerr << "Error: Spawn knot (" << spawn.x << ", " << spawn.y << ") not found in maze" << std::endl;
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
            //std::cout << line << std::endl;
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
    // Найти ID узла spawn
    int spawn_id = find_knot_id(spawn, knots);
    if (spawn_id == -1) {
        std::cerr << "Error: Spawn knot not found in knots list" << std::endl;
        return 1;
    }

    std::ifstream domain_file("DomainDatabase.json");
    json domain_data;
    domain_file >> domain_data;

    // Извлечение тегов из DomainDatabase
    std::vector<std::string> tag_list;
    for (const auto& obj : domain_data["objects"]) {
        if (obj["type"] == "knot" && obj.contains("properties") && obj["properties"].contains("tags")) {
            for (const auto& tag : obj["properties"]["tags"]) {
                tag_list.push_back(tag.get<std::string>());
            }
        }
    }

    // Желаемые дистанции
    std::map<std::string, int> desired_distances = {
        {"spawn", 0},
        {"quest_giver_heron", 5},
        {"quest_giver_cobra", 10},
        {"heron_battle", 8},
        {"cobra_hunters", 12}
    };

    // Добавление желаемых расстояний из constraints
    for (const auto& constraint : domain_data["constraints"]) {
        if (constraint.contains("desired")) {
            std::string tag = constraint["tag2"].get<std::string>();
            int desired = constraint["desired"].get<int>();
            desired_distances[tag] = desired;
        }
    }

    // Минимальные расстояния
    std::vector<std::tuple<std::string, std::string, int>> min_distance_constraints;
    for (const auto& constraint : domain_data["constraints"]) {
        if (constraint["type"] == "distance" && constraint.contains("min")) {
            min_distance_constraints.emplace_back(
                constraint["tag1"].get<std::string>(),
                constraint["tag2"].get<std::string>(),
                constraint["min"].get<int>()
            );
        }
    }

    // Запуск генетического алгоритма
    QuestPlacement quest_placement(knots, adjacencyList, tag_list, desired_distances, min_distance_constraints, spawn_id);
    auto tag_assignments = quest_placement.run_genetic_algorithm(rng);

    // Сохранение назначений
    {
        std::ofstream out("quest_assignments.txt");
        for (const auto& [tag, node_id] : tag_assignments) {
            out << "Tag: " << tag << ", Node: " << node_id << " (" << knots[node_id].x << ", " << knots[node_id].y << ")\n";
        }
    }

    return 0;
}