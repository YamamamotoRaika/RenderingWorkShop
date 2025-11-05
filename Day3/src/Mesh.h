//
// Created by nakat on 2025/10/30.
//

#ifndef MESH_H
#define MESH_H



#include "Eigen/Dense"
#include "Ray.h"
#include  "Triangle.h"

#include "Material.h"

class Mesh {
    std::vector<Triangle> triangles;
private:
   mutable int triId=-1;


public:

    Mesh() = default;
    Mesh(const std::vector<Triangle> &triangles);

    bool hit(const Ray &ray, RayHit &hit) const;
    Eigen::Vector3d getKd() const;

    bool isLight() const;
    Material getMaterial() const;
    void getNormal() const;
    Eigen::Vector3d getEmission() const;

    void setTriId(int i) const;

};



#endif //MESH_H
