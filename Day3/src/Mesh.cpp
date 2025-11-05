//
// Created by nakat on 2025/10/30.
//
#include "Ray.h"
#include "Mesh.h"
#include "Triangle.h"

void Mesh::setTriId(int i) const{
    triId = i;
};

Mesh::Mesh(const std::vector<Triangle> &tri) {
    triangles=tri;
}

bool Mesh::hit(const Ray &ray, RayHit &hit) const{
    bool isHit = false;
    hit.t = DBL_MAX;

    for (int i=0;i<triangles.size();i++) {
        Triangle tri=triangles[i];
        RayHit tempHit;
        if (tri.hit(ray, tempHit) && tempHit.t < hit.t) {
            hit = tempHit;
            isHit = true;
            setTriId(i);

        }
    }

    return isHit;
}
Eigen::Vector3d Mesh::getKd() const {
    Triangle T=triangles[triId];
    return T.getMaterial().kd * T.getMaterial().color;
}
bool Mesh::isLight() const {
    Triangle T=triangles[triId];
    return T.getMaterial().emission >0.0;
}

Eigen::Vector3d Mesh::getEmission() const{
    Triangle T= triangles[triId];
    return T.getMaterial().emission * T.getMaterial().color;
}
Material Mesh::getMaterial() const{
    Triangle T=triangles[triId];

    return  T.getMaterial();
}
