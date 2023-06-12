//
// Created by kango on 2023/03/22.
//

#include "Sphere.h"

Sphere::Sphere(double radius, Eigen::Vector3d center) : radius(radius), center(std::move(center)) {}

bool Sphere::hit(const Ray &ray, RayHit &hit) const {
    const auto b = (ray.org - center).dot(ray.dir);
    const auto c = (ray.org - center).squaredNorm() - radius * radius;
    const auto discriminant = b * b - c;

    if(discriminant < 0.0) return false;

    const Eigen::Array2d distances{(-b - sqrt(discriminant)),
                                   (-b + sqrt(discriminant))
    };

    if((distances < 1e-6).all()) return false;

    hit.t = distances[0] > 1e-6 ? distances[0] : distances[1];
    hit.point = ray.at(hit.t);
    const auto outwardNormal = (hit.point - center).normalized();

    if(outwardNormal.dot(ray.dir) < 0.0) {
        hit.normal = outwardNormal;
    } else {
        hit.normal = -outwardNormal;
    }
    return true;
}
