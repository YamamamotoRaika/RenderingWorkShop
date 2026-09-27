//
// Created by nakat on 2025/11/18.
//

#ifndef MULTILAYERBSSRDF_H
#define MULTILAYERBSSRDF_H
// MultiLayerBSSRDF.h
#pragma once
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <iostream>

struct SSSLayer {
    // 入力
    Eigen::Vector3d sigma_a;   // 吸収 [1/mm]
    Eigen::Vector3d sigma_s;   // 散乱 [1/mm]
    double g   = 0.8;          // 異方性
    double eta = 1.3;          // この層の屈折率
    double thickness = 1.0;    // 厚み [mm] (最下層は大にして半無限でもOK)

    // 事前計算（setDerived() が埋める）
    Eigen::Vector3d sigmap_s;  // σ_s' = (1-g)σ_s
    Eigen::Vector3d sigmap_t;  // σ_t' = σ_a + σ_s'
    Eigen::Vector3d D;         // 拡散係数 D = 1/(3σ_t')
    Eigen::Vector3d sig_tr;    // σ_tr = sqrt(σ_a / D)
    double Fdr = 0.0;          // 拡散フレネル（層内→層外）
    Eigen::Vector3d A;         // A = (1+Fdr)/(1-Fdr) を RGB で（Fdrは波長独立でもOK）

    void setDerived() {
        sigmap_s = (1.0 - g) * sigma_s;
        sigmap_t = sigma_a + sigmap_s;
        D        = (sigmap_t.array() * 3.0).cwiseInverse().matrix();
        // σ_tr = sqrt(σ_a / D)
        sig_tr   = (sigma_a.cwiseQuotient(D)).cwiseSqrt();
        // Fdr (η) の近似（波長独立）
        const double e = eta;
        Fdr = -1.4399/(e*e) + 0.7099/e + 0.6681 + 0.0636*e;

        //std::cout<<"Fdr = "<<Fdr<<std::endl;
        A = Eigen::Vector3d::Constant( (1.0+Fdr)/(1.0-Fdr) );

    }
    void setThickness(double l) {
        thickness=l;
    }
};

struct MLProfileSettings {
    int imagesPerSide = 8;   // 多重像の打ち切り（±m）
    bool includeSingle = false; // 単一散乱を後で足すなら true
};
Eigen::Vector3d Rd_multipole_singleLayer(
    const std::vector<SSSLayer>& layers,
    double r,
    const MLProfileSettings& opt);
Eigen::Vector3d Rd_multipole_multilayer(
    const std::vector<SSSLayer>& layers, double r, const MLProfileSettings& opt);

#endif //MULTILAYERBSSRDF_H
