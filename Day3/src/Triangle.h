//
// Created by nakat on 2025/07/09.
//

#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Ray.h"
#include "Material.h"
#include "AABB.h"

class Triangle {
    public:
    Eigen::Vector3d v0, v1, v2;
    Eigen::Vector2d uv0, uv1, uv2;
    Eigen::Vector3d n;
    Material material;

    Triangle()=default;
    Triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c);
    Triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c, const Material &material);

    Triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c,const Eigen::Vector3d &n);
    Triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c,
    const Eigen::Vector2d &d, const Eigen::Vector2d &e, const Eigen::Vector2d &f,const Eigen::Vector3d &n);

    bool hit(const Ray &ray, RayHit &hit) const;
    Material getMaterial() const;

    bool isLight() const;
    Eigen::Vector3d getKd() const;
    Eigen::Vector3d getEmission() const;
    Eigen::Vector2d getUV(const Eigen::Vector3d point) const;



    //VBH
    // 三角形の重心
    Eigen::Vector3d centroid() const;
    // 三角形のAABB
    AABB computeAABB() const;


    double area() const;
};



#endif //TRIANGLE_H
