//
// Created by nakat on 2025/11/06.
//
#include "AABB.h"

AABB::AABB()
    : min(Eigen::Vector3d::Constant(DBL_MAX)), max(Eigen::Vector3d::Constant(-DBL_MAX)) {}

AABB::AABB(const Eigen::Vector3d& min_, const Eigen::Vector3d& max_)
    : min(min_), max(max_) {}

bool AABB::hit(const Ray& ray, double tmin, double tmax) const {
    double t0 = tmin;
    double t1 = tmax;

    for (int a = 0; a < 3; ++a) {
        const double o = ray.org[a];
        const double d = ray.dir[a];

        // 平行（dir==0近傍）
        if (std::abs(d) < 1e-12) {
            // その軸のスラブ内にいなければ交差なし
            if (o < min[a] || o > max[a]) return false;
            continue;
        }

        double invD = 1.0 / d;
        double ta = (min[a] - o) * invD;
        double tb = (max[a] - o) * invD;
        if (ta > tb) std::swap(ta, tb);

        t0 = std::max(t0, ta);
        t1 = std::min(t1, tb);

        // ここは <= が安全（厚み0でも、t0==t1 で当たりとして扱いたいなら < にする）
        if (t1 < t0) return false;
    }
    return true;
}


AABB AABB::merge(const AABB& box1, const AABB& box2) {
    Eigen::Vector3d small = box1.min.cwiseMin(box2.min);
    Eigen::Vector3d big   = box1.max.cwiseMax(box2.max);
    return AABB(small, big);
}
