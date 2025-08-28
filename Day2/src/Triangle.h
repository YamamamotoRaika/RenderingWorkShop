//
// Created by nakat on 2025/07/09.
//

#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Ray.h"

class Triangle {
    public:
    Eigen::Vector3d v0, v1, v2;
    Eigen::Vector3d n;

    Triangle()=default;
    Triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c);
    Triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c,const Eigen::Vector3d &n);

    bool hit(const Ray &ray, RayHit &hit) const;
};



#endif //TRIANGLE_H
