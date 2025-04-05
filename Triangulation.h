#pragma once
//#ifndef TRIANGULATION_H
//#define TRIANGULATION_H

#include <vector>
#include <iostream>
#include <set>
#include "Knot.h"

struct Triangle {
    Knot knot_1, knot_2, knot_3;

    bool operator == (const Triangle& another) const;
};


using Edge = std::pair<Knot, Knot>;
bool operator<(const Edge& a, const Edge& b);
bool operator==(const Edge& a, const Edge& b);


double determinant(const Knot& a, const Knot& b, const Knot& c);
bool isKnotInCircumcircle(const Triangle& triangle, const Knot& knot);
std::vector<Edge> getEdgesFromTriangles(const std::vector<Triangle>& triangles);
std::vector<Triangle> delaunayTriangulation(const std::vector<Knot>& knots);
std::vector<Edge> getUniqueEdges(const std::vector<Triangle>& triangles);

//#endif // TRIANGULATION_H