//
// Created by kango on 2023/04/03.
//

#ifndef DAY_3_BODY_H
#define DAY_3_BODY_H


#include "Sphere.h"
#include "Material.h"
#include "Box.h"
#include "Triangle.h"

enum class ShapeType {
    Sphere,
    Box,
    Triangle
};

struct Body {
    ShapeType type;
    Sphere sphere;
    Box box;
    Triangle triangle;
    Material material;


    Body(Sphere sphere, Material mat) : type(ShapeType::Sphere), sphere(std::move(sphere)), material(std::move(mat)) {}
    Body(Box box, Material mat) : type(ShapeType::Box), box(std::move(box)), material(std::move(mat)) {}

    Body(Triangle triangle,Material mat) :type(ShapeType::Triangle), triangle(std::move(triangle)), material(std::move(mat)) {}

    bool hit(const Ray &ray, RayHit &hit) const {
        bool res = false;
        switch (type) {
            case ShapeType::Sphere:
                res = sphere.hit(ray, hit);
                break;
            case ShapeType::Box:
                res = box.hit(ray, hit);
                break;
            case ShapeType::Triangle:
                res = triangle.hit(ray, hit);
                break;
        }
        return res;
    }

    Eigen::Vector3d getEmission() const {
        return material.emission * material.color;
    }
    bool show() const {
        if (type == ShapeType::Triangle) {
            std::cout<< triangle.v0<<std::endl;
            std::cout<< triangle.v1<<std::endl;
            std::cout<< triangle.v2<<std::endl;

        }
        return true;
    }

    Eigen::Vector3d getKd() const {
        return material.kd * material.color;
    }

    Eigen::Vector3d getNormal(const Eigen::Vector3d &p) const;

    bool isLight() const {
        return material.emission > 0.0;
    }
};

inline Eigen::Vector3d Body::getNormal(const Eigen::Vector3d &p) const {
    switch (type) {
        case ShapeType::Sphere:
            return (p - sphere.center).normalized();
        case ShapeType::Box:
            return Eigen::Vector3d(0,0,1);
            // box.hit() で既に法線入れてる

    }
    return Eigen::Vector3d(0,0,1);
}


#endif //DAY_3_BODY_H
