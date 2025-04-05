#pragma once
#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>

struct Knot {
    float x, y;  // Изменено на float для точности
    Knot() : x(-1.0f), y(-1.0f) {}
    Knot(float _x, float _y) : x(_x), y(_y) {}

    bool operator==(const Knot& another) const;
    bool operator<(const Knot& another) const;
};

void insertKnot(std::vector<std::vector<Knot>>& grid, const Knot& knot, float cellsize);
bool isKnotValid(const std::vector<std::vector<Knot>>& grid, int g_width, int g_height,
    const Knot& knot, float radius, float width, float height, float cellsize);
float distanceBetweenPoints(float x1, float y1, float x2, float y2);
std::vector<Knot> poissonDiskSampling(float width, float height, double radius, int limit);