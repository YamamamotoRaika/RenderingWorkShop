//
// Created by nakat on 2025/11/06.
//

#ifndef AABB_H
#define AABB_H

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include "Ray.h"

class AABB {
public:
    Eigen::Vector3d min;
    Eigen::Vector3d max;

    AABB(); // デフォルトコンストラクタ
    AABB(const Eigen::Vector3d& min, const Eigen::Vector3d& max);

    // レイとの交差判定
    bool hit(const Ray& ray, double t_min = 0.001, double t_max = DBL_MAX) const;

    // マージ関数（2つのAABBを囲うAABBを返す）
    static AABB merge(const AABB& box1, const AABB& box2);
};

#endif //AABB_H
