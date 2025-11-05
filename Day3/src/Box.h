//
// Created by nakat on 2025/07/08.
//

#ifndef BOX_H
#define BOX_H

#include <vector>
#include "Triangle.h"
#include "Ray.h"

class Box {
public:
    std::vector<Triangle> triangles;

    Box() = default;

    Box(const Eigen::Vector3d &size, const Eigen::Vector3d &center);

    bool hit(const Ray &ray, RayHit &hit) const;
};

#endif



