#include <iostream>
#include <fstream>
#include <vector>
#include "Knot.h"
#include "Triangulation.h"

int main() {
    int n = 100;                // Number of knots
    float width = 100.0f;       // Maximum x-coordinate
    float height = 100.0f;      // Maximum y-coordinate
    float radius = 10;          // Minimum distance between knots
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

    std::vector<Edge> polygon = getUniqueEdges(triangles);

    {
        std::ofstream out("polygon.txt");
        for (const auto& edge : polygon) {
            out << edge.first.x << " " << edge.first.y << " "
                << edge.second.x << " " << edge.second.y << "\n";
        }
    }

    return 0;
}