//
// Created by nakat on 2025/07/08.
//
#include "Box.h"

#include "Box.h"



Box::Box(const Eigen::Vector3d &size, const Eigen::Vector3d &center) {
    Eigen::Vector3d half = size * 0.5;
    Eigen::Vector3d min = center - half;
    Eigen::Vector3d max = center + half;

    Eigen::Vector3d p000 = {min.x(), min.y(), min.z()};
    Eigen::Vector3d p001 = {min.x(), min.y(), max.z()};
    Eigen::Vector3d p010 = {min.x(), max.y(), min.z()};
    Eigen::Vector3d p011 = {min.x(), max.y(), max.z()};
    Eigen::Vector3d p100 = {max.x(), min.y(), min.z()};
    Eigen::Vector3d p101 = {max.x(), min.y(), max.z()};
    Eigen::Vector3d p110 = {max.x(), max.y(), min.z()};
    Eigen::Vector3d p111 = {max.x(), max.y(), max.z()};

    // 各面ごとに三角形を追加（12個）
    // -Z面
    triangles.emplace_back(p000, p010, p100);
    triangles.emplace_back(p100, p010, p110);
    // +Z面
    triangles.emplace_back(p001, p101, p011);
    triangles.emplace_back(p101, p111, p011);
    // -X面
    triangles.emplace_back(p000, p001, p010);
    triangles.emplace_back(p001, p011, p010);
    // +X面
    triangles.emplace_back(p100, p110, p101);
    triangles.emplace_back(p101, p110, p111);
    // -Y面
    triangles.emplace_back(p000, p100, p001);
    triangles.emplace_back(p001, p100, p101);
    // +Y面
    triangles.emplace_back(p010, p011, p110);
    triangles.emplace_back(p011, p111, p110);

}

bool Box::hit(const Ray &ray, RayHit &hit) const {
    bool isHit = false;
    hit.t = DBL_MAX;

    for (const auto& tri : triangles) {
        RayHit tempHit;
        if (tri.hit(ray, tempHit) && tempHit.t < hit.t) {
            hit = tempHit;
            isHit = true;
        }
    }

    return isHit;
}
