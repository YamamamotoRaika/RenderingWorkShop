//
// Created by nakat on 2025/11/18.
//
// MultiLayerBSSRDF.cpp
#include "MultiLayerBSSRDF.h"

#include <iostream>

// 半無限面に対する2点（実像/仮想像）からの拡散グリーン関数
// SSSLayer から全てのパラメータを取る形にする
static inline Eigen::Vector3d dipoleContribution(
    const SSSLayer& L,
    double r, double zr, double zv)
{
    const Eigen::Vector3d& sig_tr   = L.sig_tr;
    const Eigen::Vector3d& sigmap_t = L.sigmap_t;
    const Eigen::Vector3d& D        = L.D;

    auto C = [&](double rr, const Eigen::Vector3d& z) -> Eigen::Vector3d {
        Eigen::Vector3d d = Eigen::Vector3d::Constant(rr * rr) + z.cwiseProduct(z); // r^2 + z^2
        Eigen::Vector3d R = d.cwiseSqrt();
        Eigen::Vector3d e = (-sig_tr.cwiseProduct(R)).array().exp().matrix();
        // Φ = (1 / (4π D)) * e^{-σ_tr R} / R
        return e.cwiseQuotient(4.0 * EIGEN_PI * D.cwiseProduct(R));
    };

    const Eigen::Vector3d zrV = Eigen::Vector3d::Constant(zr);
    const Eigen::Vector3d zvV = Eigen::Vector3d::Constant(zv);

    Eigen::Vector3d phi_r = C(r, zrV);
    Eigen::Vector3d phi_v = C(r, zvV);

    // α' = σ_s' / σ_t'
    Eigen::Vector3d alphap = L.sigmap_s.cwiseQuotient(L.sigmap_t);

    // Rd ≈ α' ( zr Φ_r - zv Φ_v )
    return alphap.cwiseProduct( zrV.cwiseProduct(phi_r) - zvV.cwiseProduct(phi_v) );
}
// 1項の寄与: R(d;z) = z*(1+σ_tr d) / (4π d^3) * e^{-σ_tr d} t=thicknes
static inline double dipole_term(double sigma_tr, double z, double t, double d) {
    // d = sqrt(r^2 + z^2) を外で渡す
    double num = (z) * (1.0 + sigma_tr * d);
    double denom = 4.0 * EIGEN_PI * d * d * d;
    double val = num * std::exp(-sigma_tr * d) / denom ;
    return val;
}


// 有限厚 t の“multipole” : 2境界での鏡像列を ±m まで打ち切り
static inline Eigen::Vector3d Rd_singleLayer_multipole(
    const SSSLayer& L, double r, double t, int mMax)
{
    // zr, zv は RGB それぞれ
    //Eigen::Vector3d zr = Eigen::Vector3d::Ones().cwiseQuotient(L.sigmap_t);
    //Eigen::Vector3d zv = zr + 4.0 * L.A.cwiseProduct(L.D);

    Eigen::Vector3d out = Eigen::Vector3d::Zero();

    for (int c = 0; c < 3; ++c) {



        double sigma_tr = L.sig_tr[c];
        double alpha_p=L.sigmap_s[c]/L.sigmap_t[c];

        double zr = 1.0/L.sigmap_t[c];
        double zv = zr+4.0*L.A[c]*L.D[c];
        double zb=2.0*L.A[c]*L.D[c];



        double l=(L.thickness>0.0) ?2*(L.thickness + 2.0 * zb) : std::numeric_limits<double>::infinity();
        double sum = 0.0;
        auto add_pair = [&](double z_real, double z_virt, double zShift) {

            double zp =  z_real + zShift;
            double dp = std::sqrt(r*r + zp*zp);
            sum += dipole_term(sigma_tr, zp, L.thickness, dp);


            double vp = (-z_virt + zShift);
            double dvp = std::sqrt(r*r + vp*vp);
            sum -= dipole_term(sigma_tr, vp, L.thickness, dvp);
        };



        // m = ±1..mMax の鏡像列
        for (int m = -mMax; m <= mMax; ++m) {
            double shift =l*double(m);

            add_pair(zr, zv,  shift);     // 中央層
            add_pair(zr, zv,  shift + 2.0 * l); // もう一方の外挿境界からの像
        }

        out[c] = sum*alpha_p;
    }

    return out;
}


// 拡散的な“透過率”（層 i を通過して層 j に届く重み）
// まずは Beer-Lambert を σ_tr で近似。薄層で十分に効きます。
static inline Eigen::Vector3d diffuseTransmit(const SSSLayer& L, double dist){
    // dist: 層の厚さ。RGB で別々に
    return (-L.sig_tr.array()*dist).exp().matrix();
}

// N層の BSSRDF（多重散乱のみ）
// layers[0] が最上層（表面）。観測・照明は最上面から、出射も最上面。
Eigen::Vector3d Rd_multipole_multilayer(
    const std::vector<SSSLayer>& layers, double r, const MLProfileSettings& opt)
{
    const int N = (int)layers.size();
    Eigen::Vector3d Rd = Eigen::Vector3d::Zero();

    // 事前計算
    std::vector<SSSLayer> L = layers;
    for (auto& li : L) li.setDerived();


    // 各層 k に対し、その層内の multipole を計算し、
    // 上にある層を “入射で通過” × “出射で通過” で重み付けして足し込む。
    for (int k=0; k<N; ++k){
        // k層の有限厚（最下層は大きく）
        const double tk = L[k].thickness;

        // 入射で上位層を通る重み
        Eigen::Vector3d Tin = Eigen::Vector3d::Ones();
        for (int i=0; i<k; ++i) Tin = Tin.cwiseProduct( diffuseTransmit(L[i], L[i].thickness) );



        // 出射で再び通る重み
        Eigen::Vector3d Tout = Tin; // 上に戻るので同じ層を通る

        // 層内マルチポール
        Eigen::Vector3d Rk = Rd_singleLayer_multipole(L[k], r, tk, opt.imagesPerSide);

        // Fresnel 透過（最上面）をまとめて掛けたい場合はここで Ft^2 を掛ける
        Rd += Tin.cwiseProduct(Rk).cwiseProduct(Tout);


    }
    // 物理的な 1/(2π) は上位側で半径積分に合わせて調整してください
    return Rd.cwiseMax(0.0);
}

