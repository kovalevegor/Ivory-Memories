#include <iostream>
#include <fstream>
#include <vector>
#include "Knot.h"
#include "Triangulation.h"
#include "MinimumSpanningTree.h"
#include "GraphDiameter.h"
#include "SpectralClustering.h"

int main() {
    int n = 100;                // Number of knots
    float width = 192.0f;       // Maximum x-coordinate
    float height = 108.0f;      // Maximum y-coordinate
    float radius = 10;          // Minimum radius between knots
    int limit = 30;             // Maximum attempts 

    std::vector<Knot> knots;
    knots.reserve(n);

    knots = poissonDiskSampling(width, height, radius, limit);
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
    ///////////////////////////// преобразование EdgeWithWeight в Edge /////////////////////////////
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

    ///////////////////////////// Построение списка смежности /////////////////////////////
    auto adjacencyList = GraphAnalysis::buildAdjacencyList(maze);

    ///////////////////////////// Поиск диаметра графа /////////////////////////////
    auto [spawn, exit] = GraphAnalysis::findGraphDiameter(adjacencyList);

    std::cout << "Spawn knot: (" << spawn.x << ", " << spawn.y << ")\n";
    std::cout << "Exit knot: (" << exit.x << ", " << exit.y << ")\n";
    {
        std::ofstream out("spawn_exit_knots.txt");
        out << spawn.x << " " << spawn.y << "\n"
            << exit.x << " " << exit.y << "\n";
    }

    ///////////////////////////// Поиск центра графа (убежища) /////////////////////////////
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

    ///////////////////////////// Кластеризация /////////////////////////////
    int k = 7;
    float sigma = 10.0f;

    try {
        SpectralClustering sc(knots, k, sigma);
        sc.cluster();
        auto clusters = sc.getClusters();

        for (size_t i = 0; i < knots.size(); i++) {
            std::cout << "Knot (" << knots[i].x << ", " << knots[i].y
                << ") -> Cluster " << clusters[i] << std::endl;
        }

        std::ofstream out("clustering.txt");
        for (size_t i = 0; i < knots.size(); i++) {
            out << "Knot (" << knots[i].x << ", " << knots[i].y
                << ") -> Cluster " << clusters[i] << "\n";
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }


    ///////////////////////////// Распределение Квестов /////////////////////////////
    

    return 0;
}