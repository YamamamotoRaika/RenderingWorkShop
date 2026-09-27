//
// Created by kango on 2023/04/03.
//

#ifndef DAY_3_BODY_H
#define DAY_3_BODY_H


#include "Sphere.h"
#include "Material.h"
#include "Box.h"
#include "Mesh.h"
#include "Triangle.h"

enum class ShapeType {
    Sphere,
    Box,
    Triangle,
    Mesh
};

struct Body {
    ShapeType type;
    Sphere   sphere;
    Box      box;
    Mesh     mesh;
    Triangle triangle;
    Material material;


    Body(Sphere sphere, Material mat)
        : type(ShapeType::Sphere), sphere(std::move(sphere)), material(std::move(mat)) {}

    Body(Box box, Material mat)
        : type(ShapeType::Box), box(std::move(box)), material(std::move(mat)) {}

    Body(Triangle triangle, Material mat)
        : type(ShapeType::Triangle), triangle(std::move(triangle)), material(std::move(mat)) {}

    // Mesh は Mesh 側にマテリアルを持っている想定なら material は使わない
    Body(Mesh mesh)
        : type(ShapeType::Mesh), mesh(std::move(mesh)), material() {}

    bool hit(const Ray &ray, RayHit &hit) const {
        switch (type) {
            case ShapeType::Sphere:   return sphere.hit(ray, hit);
            case ShapeType::Box:      return box.hit(ray, hit);
            case ShapeType::Triangle: return triangle.hit(ray, hit);
            case ShapeType::Mesh:     return mesh.hit(ray, hit); // or hitBVH
        }
        return false;
    }
    bool Body::hit(const Ray& ray, RayHit& hit, int ignoreTriId) const {
        switch (type) {
            case ShapeType::Sphere:
                return sphere.hit(ray, hit);
            case ShapeType::Box:
                return box.hit(ray, hit);
            case ShapeType::Triangle:
                return triangle.hit(ray, hit);
            case ShapeType::Mesh:
                return mesh.hit(ray, hit, ignoreTriId);
        }
    }

    Eigen::Vector3d getEmission() const {
        // Mesh のときだけ mesh から取る
        if (type == ShapeType::Mesh) {
            return mesh.getEmission();
        }else if (type==ShapeType::Triangle) {

            return triangle.getEmission();
        }
        // それ以外は material から
        return material.emission * material.color;
    }

    Eigen::Vector3d getKd() const {
        if (type == ShapeType::Mesh) {
            return mesh.getKd();
        }else if (type==ShapeType::Triangle) {

            return triangle.getKd();
        }
        return material.kd * material.color;
    }


    Material getMaterial() const {
        if (type == ShapeType::Mesh) {
            return mesh.getMaterial();
        } else if (type==ShapeType::Triangle) {

            return triangle.getMaterial();
        }
        return material;
    }

    bool isLight() const {
        if (type == ShapeType::Mesh) {
            return mesh.isLight();
        }else if (type==ShapeType::Triangle) {

            return triangle.isLight();
        }
        return material.emission > 0.0;
    }
    int getId() const {
        if (type == ShapeType::Mesh) {
            //std::cout<<mesh.getTriId()<<std::endl;
            return mesh.getTriId();
        }
        return 0;
    }
    bool isBeard() const {
        if (type == ShapeType::Mesh) {

            return true;
        }
        return false;
    }

    Triangle getTriangle() const {

            if (type == ShapeType::Mesh) {
                //std::cout<<mesh.getTriId()<<std::endl;
                return mesh.triangles[getId()];
            }
            return Triangle();

    }

    bool show() const {
        if (type == ShapeType::Triangle) {
            std::cout << triangle.v0 << std::endl;
            std::cout << triangle.v1 << std::endl;
            std::cout << triangle.v2 << std::endl;
        } else if (type == ShapeType::Mesh) {
            std::cout << "I am mesh" << std::endl;
        }
        return true;
    }


    Eigen::Vector3d getNormal(const Eigen::Vector3d &p) const;

};

inline Eigen::Vector3d Body::getNormal(const Eigen::Vector3d &p) const {
    switch (type) {
        case ShapeType::Sphere:
            return (p - sphere.center).normalized();

        case ShapeType::Box:
            // TODO: 本当は Box の法線をちゃんと返す
            return Eigen::Vector3d(0, 0, 1);

        case ShapeType::Triangle:
            // 三角形の法線を返す（Triangle にメソッドがあるならそれを使う）
            return triangle.n; // or (v1 - v0).cross(v2 - v0).normalized();

        case ShapeType::Mesh:
            // 当面は適当な値でもいいけど、本当は mesh から三角形の法線を取る
            return Eigen::Vector3d(0, 0, 1);
    }
    return Eigen::Vector3d(0, 0, 1);
}



#endif //DAY_3_BODY_H
