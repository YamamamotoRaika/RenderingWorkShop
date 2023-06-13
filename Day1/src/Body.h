//
// Created by kango on 2023/04/03.
//

#ifndef DAY_3_BODY_H
#define DAY_3_BODY_H

#include <utility>

#include "Sphere.h"
#include "Material.h"

struct Body {
    Sphere sphere;

    Material material;

    Body(Sphere sphere, Material material) : sphere(std::move(sphere)), material(std::move(material)) {}

    bool hit(const Ray &ray, RayHit &hit) const {
        return sphere.hit(ray, hit);
    }

    Eigen::Vector3d getNormal(const Eigen::Vector3d &p) const {
        return  (p - sphere.center).normalized();
    }
};

#endif //DAY_3_BODY_H
