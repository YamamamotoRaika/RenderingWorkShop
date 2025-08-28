//
// Created by kango on 2023/04/14.
//

#ifndef DAY_3_MATERIAL_H
#define DAY_3_MATERIAL_H

#include <utility>

#include "Image.h"
#include <array>  // 忘れずに

class Material {
public:
    Color color;
    double kd;
    double emission;

    // BSSRDF用
    double eta = 1.3;
    double scale=100;
    std::array<Eigen::Vector3d, 2> sigmas;  // [0]=sigma_a, [1]=sigma_s
    bool isSubsurface = false;

    // 拡散のみ
    Material(Color color, const double &kd, const double &emission = 0.0)
        : color(std::move(color)), kd(kd), emission(emission) {}

    // BSSRDF対応コンストラクタ（型が異なるので曖昧性なし）
    Material(Color color, const double &kd, const double &emission,
             const double &eta, const Eigen::Vector3d &sigma_a,
             const Eigen::Vector3d &sigma_s, bool isSubsurface = true)
        : color(std::move(color)), kd(kd), emission(emission),
          eta(eta), sigmas{sigma_a, sigma_s}, isSubsurface(isSubsurface) {}

    // デフォルト
    Material() : color(Color(0, 0, 0)), kd(0.0), emission(0.0) {}

    Eigen::Vector3d getSigmaA() const { return sigmas[0]; }
    Eigen::Vector3d getSigmaS() const { return sigmas[1]; }
};


#endif //DAY_3_MATERIAL_H
