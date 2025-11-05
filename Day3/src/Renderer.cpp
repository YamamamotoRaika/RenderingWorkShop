//
// Created by kango on 2023/04/03.
//

#include "Renderer.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <numeric>
#include <iomanip>
#include <numbers>
#include <cmath>
#include <limits>

#ifndef TWO_PI
#define TWO_PI (double)(2.0 * EIGEN_PI)
#endif


Renderer::Renderer(const std::vector<Body> &bodies, Camera camera, Color bgColor)
        : bodies(bodies), camera(std::move(camera)), bgColor(std::move(bgColor)), engine(0), dist(0, 1) {
}

/// 乱数生成
double Renderer::rand() const {
    return dist(engine);
}

/**
 * \b シーン内に存在するBodyのうちレイにhitするものを探す
 * @param ray レイ
 * @param hit hitした物体の情報を格納するRayHit構造体
 * @return 何かしらのBodyにhitしたかどうかの真偽値
 */
bool Renderer::hitScene(const Ray &ray, RayHit &hit) const {
    /// hitするBodyのうち最小距離のものを探す
    hit.t = DBL_MAX;
    hit.idx = -1;
    for(int i = 0; i < bodies.size();++ i) {
        RayHit _hit;
        if(bodies[i].hit(ray, _hit) && _hit.t < hit.t) {
            hit.t = _hit.t;
            hit.idx = i;
            hit.point = _hit.point;
            hit.normal = _hit.normal;
        }
    }

    return hit.idx != -1;
}

Image Renderer::render() const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());
    /// フィルム上のピクセル全てに向けてレイを飛ばす
    for(int p_y = 0; p_y < image.height; p_y++) {
        for(int p_x = 0; p_x < image.width; p_x++) {
            const int p_idx = p_y * image.width + p_x;
            Color color;
            Ray ray; RayHit hit;
            camera.filmView(p_x, p_y, ray);

            /// レイを飛ばし、Bodyに当たったらその色を格納する\n
            /// 当たらなければ、背景色を返す
            //color = hitScene(ray, hit) ? bodies[hit.idx].material.color : bgColor;
            if (hitScene(ray,hit)) {
                color=bodies[hit.idx].getMaterial().color;
                //std::cout<<hit.idx<<std::endl<<" "<<bodies[hit.idx].show()<<" "<<std::endl<<std::endl;
            }else {
                color=bgColor;
            }


            image.pixels[p_idx] = color;
        }
    }

    return image;
}
double A(double eta) {
    double Fdr = -1.440 / (eta * eta) + 0.710 / eta + 0.668 + 0.0636 * eta;
    return (1 + Fdr) / (1 - Fdr);
}
static inline double FresnelDiffuseReflectance(double eta) {
    // Jensen 近似
    double invEta = 1.0 / eta, invEta2 = invEta * invEta;
    return -1.4399 * invEta2 + 0.7099 * invEta + 0.6681 + 0.0636 * eta;
}

inline double FresnelDielectric_R(double eta_i, double eta_t, double cos_i)
{
    // cos_i は [0,1] で、"外向き法線 nb" に対する入射余弦
    cos_i = std::clamp(cos_i, 0.0, 1.0);

    // Snell: sin_t^2 = (eta_i/eta_t)^2 * (1 - cos_i^2)
    double eta = eta_i / eta_t;
    double sin_t2 = eta*eta * std::max(0.0, 1.0 - cos_i*cos_i);

    // 全反射
    if (sin_t2 > 1.0) return 1.0;

    double cos_t = std::sqrt(std::max(0.0, 1.0 - sin_t2));

    double Rs = (eta_i * cos_i - eta_t * cos_t) / (eta_i * cos_i + eta_t * cos_t);
    double Rp = (eta_i * cos_t - eta_t * cos_i) / (eta_i * cos_t + eta_t * cos_i);
    return 0.5 * (Rs*Rs + Rp*Rp);
}



// 実装ファイル (Renderer.cpp など)
double  Renderer::evaluateBSSRDFScalarAtDistance(double d, const Material& material) {
    const double eta   = material.eta;
    const double scale = material.scale;
    const Eigen::Vector3d sigma_s = material.sigmas[1] * scale;
    const Eigen::Vector3d sigma_a = material.sigmas[0] * scale;
    const Eigen::Vector3d sigma_t = sigma_s + sigma_a;
    const Eigen::Vector3d sigma_tr = (sigma_a.cwiseProduct(3.0 * sigma_t)).cwiseSqrt();

    Eigen::Vector3d Rd(0,0,0);
    for (int i = 0; i < 3; ++i) {
        const double alpha_prime = sigma_s[i] / std::max(1e-16, sigma_t[i]);
        const double zr = 1.0 / std::max(1e-16, sigma_t[i]);
        const double zv = zr + 4.0 * A(eta) / (3.0 * std::max(1e-16, sigma_t[i]));
        const double dr = std::sqrt(d * d + zr * zr);
        const double dv = std::sqrt(d * d + zv * zv);

        const double phi_r = (zr * (sigma_tr[i] + 1.0 / std::max(1e-16, dr))) * std::exp(-sigma_tr[i] * dr) / std::max(1e-16, dr * dr);
        const double phi_v = (zv * (sigma_tr[i] + 1.0 / std::max(1e-16, dv))) * std::exp(-sigma_tr[i] * dv) / std::max(1e-16, dv * dv);

        Rd[i] = alpha_prime / (4.0 * EIGEN_PI) * (phi_r + phi_v);
    }
    // サンプリングは輝度でなくてもOK。簡潔に RGB 平均を採用（必要なら ITU 加重に差し替え可）
    return (Rd[0] + Rd[1] + Rd[2]) / 3.0;
}
double Renderer::evaluateBSSRDFScalarAtDistanceSS(double r, const Material& material) {
    const double scale = material.scale;
    // あなたの設計では sigmas[1] を σs′（reduced）として扱っている前提
    const Eigen::Vector3d sigma_a   = material.sigmas[0] * scale;
    const Eigen::Vector3d sigma_s_p = material.sigmas[1] * scale;

    // RGB をスカラー（平均）にまとめる
    double sum = 0.0;
    for (int c = 0; c < 3; ++c) {
        const double sigma_t_p = std::max(1e-16, sigma_a[c] + sigma_s_p[c]); // σt′
        const double alpha_p   = std::max(0.0, std::min(1.0, sigma_s_p[c] / sigma_t_p));
        const double z_e       = 1.0 / sigma_t_p;                              // 有効深さ（近似）
        const double dr        = std::sqrt(r*r + z_e*z_e);
        const double att       = std::exp(-2.0 * sigma_t_p * dr);              // 入出のビール減衰
        const double Fdr       = std::clamp(FresnelDiffuseReflectance(material.eta), 0.0, 0.999);
        const double Ft        = 1.0 - Fdr;                                    // 透過（入）近似
        // 係数：Ft^2（入・出）、α′、幾何(1/4π)、距離(1/dr^2)
        const double C         = (Ft * Ft) * alpha_p / (4.0 * EIGEN_PI);
        sum += C * att / std::max(1e-16, dr * dr);
    }
    return sum / 3.0; // スカラー（平均）
}

// ---- Fresnel の拡散反射率（Jensen近似）----
static inline double Fdr_from_eta(double eta) {
    // eta>=1 の近似（eta<1 でも誤差は小）
    double invEta = 1.0 / eta;
    double invEta2 = invEta * invEta;
    return -1.4399*invEta2 + 0.7099*invEta + 0.6681 + 0.0636*eta;
}

// 1項の寄与: R(d;z) = z*(1+σ_tr d) / (4π d^3) * e^{-σ_tr d} t=thicknes
static inline double dipole_term(double sigma_tr, double z, double t, double d) {
    // d = sqrt(r^2 + z^2) を外で渡す
    double num = (z) * (1.0 + sigma_tr * d);
    double denom = 4.0 * EIGEN_PI * d * d * d;
    double val = num * std::exp(-sigma_tr * d) / denom ;
    return val;
}

static inline double _dipole_term(double sigma_tr, double z, double t, double d) {
    // d = sqrt(r^2 + z^2) を外で渡す
    double num = (t-z) * (1.0 + sigma_tr * d);
    double denom = 4.0 * EIGEN_PI * d * d * d;
    double val = num * std::exp(-sigma_tr * d) / denom ;
    return std::max(0.0, val);
}

// r_max: 比較に使う最大半径（例: 0.02 m = 20 mm）
// zb, T は理論計算内で求めた値
int chooseImageOrder(double r_max, double L) {
    // r_max から ±mL まで像を含める目安：両側合わせて ~ (r_max/L)*2 + マージン
    const int m = int(std::ceil( (r_max / L) * 2.0 )) + 3;
    return std::max(5, std::min(m, 100));  // 過剰振動防止に上限も
}


// r だけに依存するラジアルプロファイル（RGBベクトル）
static inline Eigen::Vector3d Rd_multipole_RGB(const Eigen::Vector3d& sigma_a,
                                               const Eigen::Vector3d& sigma_sp, // σ_s'
                                               double A, double r,
                                               double rMax,
                                               double thickness)
{
    Eigen::Vector3d out = Eigen::Vector3d::Zero();

    for (int c = 0; c < 3; ++c) {
        double sa = std::max(1e-9, sigma_a[c]);
        double ss = std::max(1e-9, sigma_sp[c]);
        double st = sa + ss;                           // σ_t'
        double alpha_p = ss / st;                      // α'
        double D = 1.0 / (3.0 * st);                   // 1/(3σ_t')
        double sigma_tr = std::sqrt(3.0 * sa * st);    // √(3σ_aσ_t')
        double zr = 1.0 / st;                          // 実源
        double zv = zr + 4.0 * A * D;                  // 仮想源（半無限）
        double zb=2*A*D;

        //std::cout<<"zr="<<zr<<",zv="<<zv<<std::endl;

        // 有限厚のときは 2L 周期で像を足す（上下境界に無限像）
        // L = 2(T + 4AD) （外挿境界までの距離）
        double L = (thickness > 0.0) ? 2*(thickness + 2.0 * zb) : std::numeric_limits<double>::infinity();

        //const int M=chooseImageOrder(rMax,L);
        const int M=3;


        double sum = 0.0;
        auto add_pair = [&](double z_real, double z_virt, double zShift) {

            double zp =  z_real + zShift;
            double dp = std::sqrt(r*r + zp*zp);

            sum += dipole_term(sigma_tr, zp, thickness, dp);



            double vp = (-z_virt + zShift);
            double dvp = std::sqrt(r*r + vp*vp);

            sum -= dipole_term(sigma_tr, vp, thickness, dvp);

        };

        if (std::isinf(L)) {
            std::cout<<"dipole"<<std::endl;
            // 半無限：通常のディポール（上下対称を2倍カウント）
            double dr = std::sqrt(r*r + zr*zr);
            double dv = std::sqrt(r*r + zv*zv);
            sum  = dipole_term(sigma_tr, zr, r, dr);
            sum += dipole_term(sigma_tr, zv, r, dv);
            sum *= 2.0; // ±z を合算（add_pair 相当）
        } else {
            // マルチポール：像の列を -M..M で切る（M=2〜3で十分）
            // 周期 2L ごとに上下へ平行移動した像を追加
            //std::cout<<"multipole,L="<< L<<"thickness=" <<thickness<<std::endl;
            for (int m = -M; m <= M; ++m) {
                double shift = L * double(m);
                add_pair(zr, zv,  shift);     // 中央層
                add_pair(zr, zv,  shift + 2.0 * L); // もう一方の外挿境界からの像
            }
        }

        out[c] = alpha_p * sum;
    }


    // clamp
    out = out.cwiseMax(0.0);
    return out;
}

double Renderer::evaluateBSSRDF_MultipoleAtDistance(const double r,
                                          const Material& material,
                                          double thickness,
                                          double r_max,
                                          int    M)
{
    // r のみで決まる近似（平面局所近似）


    // マテリアル：あなたの約束に合わせる
    // sigmas[0]=σ_a, sigmas[1]=σ_s' （すでに reduced scattering で持っている前提）
    double g=0;
    Eigen::Vector3d sigma_a = material.sigmas[0] * material.scale;
    Eigen::Vector3d sigma_sp = (1.0 - g) * material.sigmas[1] * material.scale;
    // Fresnel の拡散反射率→A
    const double Fdr = std::clamp(Fdr_from_eta(material.eta), 0.0, 0.99);
    const double a   = A(material.eta);


    Eigen::Vector3d Rd = Rd_multipole_RGB(sigma_a, sigma_sp, a, r,r_max, thickness);

    // 入出の平均フレネル透過（角度平均）と 1/π を掛けて BSSRDF へ
    const double Ft_bar = 1.0 - Fdr;          // 角度平均近似
    const double scaleF = 1 / EIGEN_PI;
    Rd *= scaleF;

    // Rd は [1/area] のラジアルプロファイル（RGB）
    return (Rd[0]+ Rd[1] + Rd[2]) /3.0;
}


// 実装
void BSSRDFRadialCDF::build(const Material& material, int N) {
    // 代表的な減衰率：RGB の sigma_tr の平均で Rmax を決める（10/σ_tr 目安）
    const double scale = material.scale;
    Eigen::Vector3d sigma_s = material.sigmas[1] * scale;
    Eigen::Vector3d sigma_a = material.sigmas[0] * scale;
    Eigen::Vector3d sigma_t = sigma_s + sigma_a;
    Eigen::Vector3d sigma_tr = (sigma_a.cwiseProduct(3.0 * sigma_t)).cwiseSqrt();
    const double sigma_tr_avg = std::max(1e-12, (sigma_tr[0] + sigma_tr[1] + sigma_tr[2]) / 3.0);

    rMax = 10.0 / sigma_tr_avg;
    //rMax = 2.0;

    rGrid.resize(N);
    cdf.resize(N);
    double sum = 0.0;
    double prev_r = 0.0;
    double prev_val = 0.0; // integrand at r=0 は 0

    for (int i = 0; i < N; ++i) {
        const double r = rMax * double(i) / double(N - 1);
        rGrid[i] = r;
        //const double Rd = Renderer::evaluateBSSRDFScalarAtDistance(r, material);
        const double Rd = Renderer::evaluateBSSRDF_MultipoleAtDistance(r, material,200,0.02,5);


        const double val = 2.0 * EIGEN_PI * r * std::max(0.0, Rd); // integrand
        if (i > 0) sum += 0.5 * (val + prev_val) * (r - prev_r);    // 台形則
        cdf[i] = sum;
        prev_r = r; prev_val = val;
    }
    Z = std::max(1e-18, sum);
    for (int i = 0; i < N; ++i) cdf[i] /= Z;

    // 尾部が残っていれば一度だけ拡張
    if (cdf.back() < 0.999) {
        rMax *= 2.0;
        build(material, N);
    }
}
// BSSRDFRadialCDF::buildSS 実装（Renderer.cpp 側に書くのが楽）
void BSSRDFRadialCDF::buildSS(const Material& material, int N) {
    // MS の build と同様。違いは Rd 呼び先だけ
    const double scale = material.scale;
    Eigen::Vector3d sigma_s = material.sigmas[1] * scale;
    Eigen::Vector3d sigma_a = material.sigmas[0] * scale;
    Eigen::Vector3d sigma_t = sigma_s + sigma_a;
    Eigen::Vector3d sigma_tr = (sigma_a.cwiseProduct(3.0 * sigma_t)).cwiseSqrt();
    const double sigma_tr_avg = std::max(1e-12, (sigma_tr[0] + sigma_tr[1] + sigma_tr[2]) / 3.0);

    rMax = 10.0 / sigma_tr_avg;

    rGrid.resize(N);
    cdf.resize(N);
    double sum = 0.0, prev_r = 0.0, prev_val = 0.0;

    for (int i = 0; i < N; ++i) {
        const double r = rMax * double(i) / double(N - 1);
        rGrid[i] = r;
        const double Rd = Renderer::evaluateBSSRDFScalarAtDistanceSS(r, material); // ★ SS
        const double val = 2.0 * EIGEN_PI * r * std::max(0.0, Rd);
        if (i > 0) sum += 0.5 * (val + prev_val) * (r - prev_r);
        cdf[i] = sum;
        prev_r = r; prev_val = val;
    }
    Z = std::max(1e-18, sum);
    for (int i = 0; i < N; ++i) cdf[i] /= Z;

    if (cdf.back() < 0.999) { rMax *= 2.0; buildSS(material, N); }
}


double BSSRDFRadialCDF::sampleR(double u) const {
    // 2分探索＋線形補間
    auto it = std::lower_bound(cdf.begin(), cdf.end(), u);
    int idx = int(it - cdf.begin());
    if (idx <= 0) return rGrid.front();
    if (idx >= int(cdf.size())) return rGrid.back();
    const double c0 = cdf[idx - 1], c1 = cdf[idx];
    const double t  = (u - c0) / std::max(1e-16, (c1 - c0));
    return rGrid[idx - 1] + t * (rGrid[idx] - rGrid[idx - 1]);
}


// Renderer.cpp
const BSSRDFRadialCDF& Renderer::getRadialTable(const Material& m) const {
    // すでにあるなら即返す
    {
        std::lock_guard<std::mutex> lk(sssCacheMtx);
        auto it = sssRadialCache.find(&m);
        if (it != sssRadialCache.end()) return *(it->second);
    }
    // なければ作る（ロック外で重い処理）
    auto tbl = std::make_shared<BSSRDFRadialCDF>();
    tbl->build(m);

    // 登録して返す
    {
        std::lock_guard<std::mutex> lk(sssCacheMtx);
        auto [it, _] = sssRadialCache.emplace(&m, std::move(tbl));
        return *(it->second);
    }
}
// Renderer.cpp に追加：MS/SS それぞれビルド（中身は build() と同じで Rd 呼び先だけ違う）
const BSSRDFRadialCDF& Renderer::getRadialTableMS(const Material& m) const {
    { std::lock_guard<std::mutex> lk(sssCacheMtx);
        auto it = sssRadialCacheMS.find(&m);
        if (it != sssRadialCacheMS.end()) return *(it->second); }
    auto tbl = std::make_shared<BSSRDFRadialCDF>();
    tbl->build(m); // ← 既存（MS: evaluateBSSRDFScalarAtDistance を使う）
    { std::lock_guard<std::mutex> lk(sssCacheMtx);
        auto [it,_] = sssRadialCacheMS.emplace(&m, std::move(tbl));
        return *(it->second); }
}

const BSSRDFRadialCDF& Renderer::getRadialTableSS(const Material& m) const {
    { std::lock_guard<std::mutex> lk(sssCacheMtx);
        auto it = sssRadialCacheSS.find(&m);
        if (it != sssRadialCacheSS.end()) return *(it->second); }
    auto tbl = std::make_shared<BSSRDFRadialCDF>();
    // ★ SS 用に buildSS を用意して呼ぶ
    tbl->buildSS(m);
    { std::lock_guard<std::mutex> lk(sssCacheMtx);
        auto [it,_] = sssRadialCacheSS.emplace(&m, std::move(tbl));
        return *(it->second); }
}


Color Renderer::evaluateBSSRDF_Multipole(const Eigen::Vector3d& xi,
                                          const Eigen::Vector3d& xo,
                                          const Material& material,
                                          double thickness,
                                          double r_max,
                                          int    M) const
{
    // r のみで決まる近似（平面局所近似）
    const double r = (xo - xi).norm();

    // マテリアル：あなたの約束に合わせる
    // sigmas[0]=σ_a, sigmas[1]=σ_s' （すでに reduced scattering で持っている前提）
    double g=0;
    Eigen::Vector3d sigma_a = material.sigmas[0] * material.scale;
    Eigen::Vector3d sigma_sp = (1.0 - g) * material.sigmas[1] * material.scale;
    // Fresnel の拡散反射率→A
    const double Fdr = std::clamp(Fdr_from_eta(material.eta), 0.0, 0.99);
    const double a   = A(material.eta);


    Eigen::Vector3d Rd = Rd_multipole_RGB(sigma_a, sigma_sp, a, r,r_max, thickness);

    // 入出の平均フレネル透過（角度平均）と 1/π を掛けて BSSRDF へ
    const double Ft_bar = 1.0 - Fdr;          // 角度平均近似
    const double scaleF = Ft_bar / EIGEN_PI;
    //Rd *= scaleF;

    // Rd は [1/area] のラジアルプロファイル（RGB）
    return Color(Rd[0], Rd[1], Rd[2]);
}




Image Renderer::directIlluminationRender(const unsigned int &samples) const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());
    /// フィルム上のピクセル全てに向けてレイを飛ばす
#pragma omp parallel for
    for(int p_y = 0; p_y < image.height; p_y++) {
        for(int p_x = 0; p_x < image.width; p_x++) {
            const int p_idx = p_y * image.width + p_x;
            Ray ray; RayHit hit;
            camera.filmView(p_x, p_y, ray);

            if(hitScene(ray, hit)) {
                //debug用
                //hit.show();
                Color reflectRadiance = Color::Zero();
                for(int i = 0; i < samples; ++i) {
                    /// 衝突点xから半球上のランダムな方向にレイを飛ばす
                    Ray _ray; RayHit _hit;
                    diffuseSample(hit.point, hit.normal, _ray);

                    /// もしBodyに当たったら,その発光量を加算する
                    if(hitScene(_ray, _hit)) {


                        reflectRadiance += bodies[hit.idx].getKd().cwiseProduct(bodies[_hit.idx].getEmission());
                    }
                }
                /// 自己発光 + 反射光
                image.pixels[p_idx] = bodies[hit.idx].getEmission() + reflectRadiance / static_cast<double>(samples);
            } else {
                image.pixels[p_idx] = bgColor;
            }

        }
    }

    return image;
}

Image Renderer::_directIlluminationRender(const unsigned int &samples) const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());
    /// フィルム上のピクセル全てに向けてレイを飛ばす
#pragma omp parallel for
    for(int p_y = 0; p_y < image.height; p_y++) {
        for(int p_x = 0; p_x < image.width; p_x++) {
            const int p_idx = p_y * image.width + p_x;
            Ray ray;
            RayHit hit;
            camera.filmView(p_x, p_y, ray);

            if (hitScene(ray, hit)) {
                if(bodies[hit.idx].isLight()) {
                    image.pixels[p_idx] = bodies[hit.idx].getEmission();
                } else {
                    Color reflectRadiance = Color::Zero();
                    for (int i = 0; i < samples; ++i) {
                        /// 衝突点hit.pointから半球上のランダムな方向にレイを飛ばす
                        Ray _ray; RayHit _hit;
                        diffuseSample(hit.point, hit.normal, _ray);

                        /// もしBodyに当たったら,その発光量を加算する
                        if (hitScene(_ray, _hit) && bodies[_hit.idx].isLight()) {
                            reflectRadiance += bodies[hit.idx].getKd().cwiseProduct(bodies[_hit.idx].getEmission());
                        }
                    }
                    /// 自己発光 + 反射光(今回、光源以外に自己発光している物体はなく、光源の場合は除外しているので自己発光の部分は梨)
                    image.pixels[p_idx] = reflectRadiance / static_cast<double>(samples);
                }
            } else {
                image.pixels[p_idx] = bgColor;
            }
        }
    }

    return image;
}

void Renderer::diffuseSample(const Eigen::Vector3d &incidentPoint, const Eigen::Vector3d &normal, Ray &out_Ray) const {
    const double phi = 2.0 * EIGEN_PI * rand();
    const double theta = asin(sqrt(rand()));

    /// normalの方向をy軸とした正規直交基底を作る
    Eigen::Vector3d u, v;
    computeLocalFrame(normal, u, v);

    const double _x = sin(theta) * cos(phi);
    const double _y = cos(theta);
    const double _z = sin(theta) * sin(phi);

    out_Ray.dir = _x * u + _y * normal + _z * v;
    out_Ray.org = incidentPoint;
}
// Fresnel（Schlick近似でOK）
inline double fresnelT_schlick(double eta, double cosTheta) {
    cosTheta = std::clamp(cosTheta, 0.0, 1.0);
    const double R0 = std::pow((1.0 - eta) / (1.0 + eta), 2.0);
    const double Fr = R0 + (1.0 - R0) * std::pow(1.0 - cosTheta, 5.0);
    return 1.0 - Fr; // 透過 Ft
}

//radiusの値って具体的になんだろう？
void Renderer::SSSSample(const Eigen::Vector3d &incidentPoint, const double &radius,double &r,const Eigen::Vector3d &normal, Ray &out_Ray) const {
    const double phi = 2.0 * EIGEN_PI * rand();
    const double theta = asin(sqrt(rand()));
    double ran=(double)rand() / (double)RAND_MAX;
    r=radius*ran;


    /// normalの方向をy軸とした正規直交基底を作る
    Eigen::Vector3d u, v;
    computeLocalFrame(normal, u, v);

    const double _x = sin(theta) * cos(phi);
    const double _y = cos(theta);
    const double _z = sin(theta) * sin(phi);

    //incident_pointから距離_r離した点を考える
    //角度は今回Φを流用する。修正筆頭ポイント
    const double d_x =r*cos(phi);
    const double d_y =r*sin(phi);

    out_Ray.org = incidentPoint+d_x*u+d_y*v;
}

// 変更前： void Renderer::_SSSSample(..., Ray &out_Ray) const
// 変更後： pdfA を返す引数を追加
void Renderer::_SSSSample(const Eigen::Vector3d &incidentPoint,
                          const Eigen::Vector3d &normal,
                          const Material& material,
                          Ray &out_Ray,
                          double &out_pdfA) const
{
    const auto& tbl = getRadialTable(material);


    const double u1 = (double)rand() / (double)RAND_MAX;
    const double u2 = rand();
    const double r  = tbl.sampleR(u1);
    const double phi = 2.0 * EIGEN_PI * u2;

    Eigen::Vector3d u, v;
    computeLocalFrame(normal, u, v);

    const double dx = r * std::cos(phi);
    const double dy = r * std::sin(phi);
    out_Ray.org = incidentPoint + dx * u + dy * v;

    // ★ 面積pdf： p_A = Rd(r) / Z
    //const double Rd = evaluateBSSRDFScalarAtDistance(r,material);
    const Color RdC = evaluateBSSRDF_Multipole(incidentPoint,out_Ray.org,material,200,0.2,5);
    const double Rd=(RdC[1]+RdC[2]+RdC[0])/3.0;
    out_pdfA = std::max(1e-18, Rd / tbl.Z);

    //out_pdfA = std::max(1e-18, Rd / tbl.Z);
}
void Renderer::_SSSSampleMixed(const Eigen::Vector3d &incidentPoint,
                               const Eigen::Vector3d &normal,
                               const Material& material,
                               Ray &out_Ray,
                               double &out_pdfA) const
{
    const auto& ms = getRadialTableMS(material);
    const auto& ss = getRadialTableSS(material);

    const double Zms = ms.Z;
    const double Zss = ss.Z;
    const double Zsum= std::max(1e-18, Zms + Zss);
    const double wMS = Zms / Zsum;
    const double wSS = Zss / Zsum;

    // 成分選択
    const double uComp = rand();
    const bool chooseMS = (uComp < wMS);
    const BSSRDFRadialCDF& tbl = chooseMS ? ms : ss;

    // r, φ サンプル
    const double u1 = rand(), u2 = rand();
    const double r  = tbl.sampleR(u1);
    const double phi= 2.0 * EIGEN_PI * u2;

    Eigen::Vector3d t,b; computeLocalFrame(normal, t, b);
    out_Ray.org = incidentPoint + (r*std::cos(phi))*t + (r*std::sin(phi))*b;

    // 混合 面積pdf
    const double Rd_ms = Renderer::evaluateBSSRDFScalarAtDistance(r, material);
    const double Rd_ss = Renderer::evaluateBSSRDFScalarAtDistanceSS(r, material);
    const double pdfA_ms = (Zms > 0.0) ? Rd_ms / Zms : 0.0;
    const double pdfA_ss = (Zss > 0.0) ? Rd_ss / Zss : 0.0;
    out_pdfA = std::max(1e-18, wMS * pdfA_ms + wSS * pdfA_ss);
}





Color Renderer::evaluateBSSRDF(const Eigen::Vector3d& xi,
                                const Eigen::Vector3d& xo,
                                const Material& material) const {
    const double eta   = material.eta;
    const double scale = material.scale;
    const Eigen::Vector3d sigma_s = material.sigmas[1] * scale;
    const Eigen::Vector3d sigma_a = material.sigmas[0] * scale;
    const Eigen::Vector3d sigma_t = sigma_s + sigma_a;
    const Eigen::Vector3d sigma_tr = (sigma_a.cwiseProduct(3.0 * sigma_t)).cwiseSqrt();

    const double d = (xo - xi).norm();  // 距離

    Eigen::Vector3d Rd;

    for (int i = 0; i < 3; ++i) {
        // 多くの論文で使われる近似（dipoleプロファイル）
        const double alpha_prime = sigma_s[i] / sigma_t[i];
        const double zr = 1.0 / sigma_t[i];
        const double zv = zr + 4.0 * A(eta) / (3.0 * sigma_t[i]);
        const double dr = std::sqrt(d * d + zr * zr);
        const double dv = std::sqrt(d * d + zv * zv);

        const double phi_r = (zr * (sigma_tr[i] + 1.0 / dr)) * std::exp(-sigma_tr[i] * dr) / (dr * dr);
        const double phi_v = (zv * (sigma_tr[i] + 1.0 / dv)) * std::exp(-sigma_tr[i] * dv) / (dv * dv);

        Rd[i] = alpha_prime / (4.0 * EIGEN_PI) * (phi_r + phi_v);
    }
   // std::cout << "Rd:"<< Rd[0]<<":"<<Rd[1]<<":"<<Rd[2]<< std::endl;
   // std::cout << std::endl;
    return Rd;  // RGBに対応した表面下拡散反射
}
Color Renderer::_evaluateBSSRDF(const Eigen::Vector3d& xi,
                               const Eigen::Vector3d& xo,
                               const Material& material) const {
    const double d = (xo - xi).norm();

    // 既存 MS（あなたの実装のままでOK）
    Color Rd_ms =evaluateBSSRDF(xi, xo, material);

    // SS を RGB で加算（単一散乱近似）
    const double scale = material.scale;
    const Eigen::Vector3d sigma_a   = material.sigmas[0] * scale;
    const Eigen::Vector3d sigma_s_p = material.sigmas[1] * scale;
    Eigen::Vector3d Rd_ss_rgb(0,0,0);
    const double Fdr = std::clamp(FresnelDiffuseReflectance(material.eta), 0.0, 0.999);
    const double Ft  = 1.0 - Fdr;
    for (int c = 0; c < 3; ++c) {
        const double sigma_t_p = std::max(1e-16, sigma_a[c] + sigma_s_p[c]);
        const double alpha_p   = std::max(0.0, std::min(1.0, sigma_s_p[c] / sigma_t_p));
        const double z_e       = 1.0 / sigma_t_p;
        const double dr        = std::sqrt(d*d + z_e*z_e);
        const double att       = std::exp(-2.0 * sigma_t_p * dr);
        const double C         = (Ft * Ft) * alpha_p / (4.0 * EIGEN_PI);
        Rd_ss_rgb[c] = C * att / std::max(1e-16, dr*dr);
    }

    return Rd_ms + Rd_ss_rgb; // MS + SS の合成
}



Image Renderer::SSSdirectIlluminationRender(const unsigned int &samples) const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());

#pragma omp parallel for
    for(int p_y = 0; p_y < image.height; p_y++) {
        for(int p_x = 0; p_x < image.width; p_x++) {
            const int p_idx = p_y * image.width + p_x;
            Ray ray;
            RayHit hit;
            camera.filmView(p_x, p_y, ray);

            if (hitScene(ray, hit)) {
                if (bodies[hit.idx].isLight()) {
                    image.pixels[p_idx] = bodies[hit.idx].getEmission();
                } else {
                    Color reflectRadiance = Color::Zero();

                    for (int i = 0; i < samples; ++i) {
                        // … ピクセル内ループ …
                        Ray _ray; RayHit _hit;

                        // ① 位置：Rd に沿って
                        double pdfA = 0.0;
                        _SSSSample(hit.point, hit.normal, bodies[hit.idx].material, _ray, pdfA); // xi = _ray.org

                        // ② 方向：既存の diffuseSample を使う（コサイン分布）
                        //    → 後で p_dir = cosθ/π で割る
                        diffuseSample(_ray.org, hit.normal, _ray);
                        const double cos_i = std::max(0.0, hit.normal.dot(_ray.dir));



                        if (hitScene(_ray, _hit) && bodies[_hit.idx].isLight()) {
                            const Color Li = bodies[_hit.idx].getEmission();

                            // ③ S = (Ft_i Ft_o / π) * Rd_rgb(d)
                            const double cos_o = std::max(0.0, hit.normal.dot(-ray.dir)); // カメラ側
                            const double Ft_i  = fresnelT_schlick(bodies[hit.idx].material.eta, cos_i);
                            const double Ft_o  = fresnelT_schlick(bodies[hit.idx].material.eta, cos_o);

                            const Color Rd_rgb = evaluateBSSRDF(_ray.org, hit.point, bodies[hit.idx].material);
                            const Color S_full = (Ft_i * Ft_o / EIGEN_PI) * Rd_rgb;

                            //std::cout << (Ft_i * Ft_o / EIGEN_PI)<<std::endl;


                            // ④ 不偏推定：
                            reflectRadiance += (Ft_i * Ft_o )*bodies[hit.idx].getKd().cwiseProduct( Li ) / (pdfA );
                        }
                    }

                    image.pixels[p_idx] = reflectRadiance / static_cast<double>(samples);
                }
            } else {
                image.pixels[p_idx] = bgColor;
            }
        }
    }

    return image;
}
Image Renderer::_SSSdirectIlluminationRender(const unsigned int &samples) const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());

#pragma omp parallel for
    for(int p_y = 0; p_y < image.height; p_y++) {
        for(int p_x = 0; p_x < image.width; p_x++) {
            const int p_idx = p_y * image.width + p_x;
            Ray ray;
            RayHit hit;
            camera.filmView(p_x, p_y, ray);

            if (hitScene(ray, hit)) {
                if (bodies[hit.idx].isLight()) {
                    image.pixels[p_idx] = bodies[hit.idx].getEmission();
                } else {
                    Color reflectRadiance = Color::Zero();

                    for (int i = 0; i < samples; ++i) {
                        Ray _ray; RayHit _hit;

                        // ★ 基準点は hit.point を渡す（未初期化の _ray.org ではない）
                        double pdfA = 0.0;
                        _SSSSample(hit.point, hit.normal, bodies[hit.idx].material, _ray, pdfA);

                        // 出射方向のサンプル（従来通り）
                        diffuseSample(_ray.org, hit.normal, _ray);

                        if (hitScene(_ray, _hit) && bodies[_hit.idx].isLight()) {
                            const Color Li = bodies[_hit.idx].getEmission();                 // 光源放射
                            const Color S  = evaluateBSSRDF(_ray.org, hit.point, bodies[hit.idx].material); // Rd評価（RGB）
                            // ★ 面積pdfで割る
                           // reflectRadiance += bodies[hit.idx].getKd().cwiseProduct(S.cwiseProduct(Li)) / pdfA;
                            reflectRadiance += bodies[hit.idx].getKd().cwiseProduct(Li) / pdfA ;
                        }
                    }

                    image.pixels[p_idx] = reflectRadiance / static_cast<double>(samples);
                }
            } else {
                image.pixels[p_idx] = bgColor;
            }
        }
    }

    return image;
}// ======================= Renderer.cpp 差し替えブロック =======================
// 出口点を実表面に投影（同一 Body 上を保証）
bool Renderer::projectToSurface(const Eigen::Vector3d& xop,
                                const Eigen::Vector3d& n_hint,
                                int bodyId,
                                Eigen::Vector3d& xo,
                                Eigen::Vector3d& n_o) const {
    const double eps = 1e-4; // シーンスケールに合わせて調整可
    for (int s = 0; s < 2; ++s) {
        const Eigen::Vector3d dir = (s == 0) ?  n_hint : -n_hint;
        Ray r; r.org = xop + dir * eps; r.dir = dir;
        RayHit h;
        if (hitScene(r, h) && h.idx == bodyId) {
            xo  = h.point;
            n_o = h.normal;
            return true;
        }
    }
    return false;
}

// 反射率（フレネル）Schlick 近似（透過係数 Ft = 1 - Fr）
static inline double fresnelReflectanceSchlick(double eta, double cosTheta) {
    // cosTheta は入射側の面法線との内積の絶対値（[0,1]）
    cosTheta = std::clamp(std::abs(cosTheta), 0.0, 1.0);
    // 入外どちらにいても概形が崩れないよう Schlick 近似で簡略化
    // R0 = ((n1-n2)/(n1+n2))^2。ここでは n1=1, n2=eta の想定。
    const double r0 = (1.0 - eta) / (1.0 + eta);
    const double R0 = r0 * r0;
    return R0 + (1.0 - R0) * std::pow(1.0 - cosTheta, 5.0);
}
static inline double fresnelTransmitSchlick(double eta, double cosTheta) {
    return 1.0 - fresnelReflectanceSchlick(eta, cosTheta); // Ft ≈ 1 - Fr
}



// ========== 乱数ユーティリティ ==========
inline double Renderer::urand() const {
    // あなたのプロジェクトの RNG に合わせてください。
    // ここでは C++ の thread_local Xorshift もどき
    thread_local uint64_t s = 88172645463393265ull;
    s ^= s << 7; s ^= s >> 9;
    return (s * 0x9E3779B97F4A7C15ull) * (1.0 / double(UINT64_MAX));
}

// ========== 幾何ユーティリティ ==========
Eigen::Vector3d Renderer::sampleIsotropicDir(double u1, double u2) const {
    const double z  = 1.0 - 2.0 * u1;                 // cosθ ∈ [-1,1]
    const double r  = std::sqrt(std::max(0.0, 1.0 - z*z));
    const double phi= 2.0 * EIGEN_PI * u2;
    return Eigen::Vector3d(r*std::cos(phi), r*std::sin(phi), z).normalized();
}

double Renderer::sampleFreeFlight(double u, double sigma_t) const {
    if (sigma_t <= 0.0) return std::numeric_limits<double>::infinity();
    return -std::log(1.0 - std::max(0.0, std::min(1.0, u))) / sigma_t;
}

bool Renderer::refract(const Eigen::Vector3d &wi, const Eigen::Vector3d &n, double n1, double n2,
                       Eigen::Vector3d &wt) const {
    // wi: 入射（外向き基準、表面に向かうときは - を付けて渡す）
    double eta = n1 / n2;
    double cosI = -n.dot(wi);
    double sin2T = eta*eta * std::max(0.0, 1.0 - cosI*cosI);
    if (sin2T > 1.0) return false; // 全反射
    double cosT = std::sqrt(std::max(0.0, 1.0 - sin2T));
    wt = eta * wi + (eta * cosI - cosT) * n;
    wt.normalize();
    return true;
}

double Renderer::fresnelDielectric_R(double n1, double n2, double cosThetaI) const {
    // 完全版：スネルから cosT を計算して Fresnel 反射率を返す
    cosThetaI = std::clamp(std::abs(cosThetaI), 0.0, 1.0);
    double etaI = n1, etaT = n2;
    bool entering = true;
    if (cosThetaI < 0.0) { entering = false; std::swap(etaI, etaT); cosThetaI = std::abs(cosThetaI); }
    double sinThetaI = std::sqrt(std::max(0.0, 1.0 - cosThetaI*cosThetaI));
    double sinThetaT = etaI / etaT * sinThetaI;
    if (sinThetaT >= 1.0) return 1.0; // 全反射
    double cosThetaT = std::sqrt(std::max(0.0, 1.0 - sinThetaT*sinThetaT));
    double r_parl = ((etaT * cosThetaI) - (etaI * cosThetaT)) / ((etaT * cosThetaI) + (etaI * cosThetaT));
    double r_perp = ((etaI * cosThetaI) - (etaT * cosThetaT)) / ((etaI * cosThetaI) + (etaT * cosThetaT));
    return 0.5 * (r_parl*r_parl + r_perp*r_perp);
}

inline Eigen::Vector3d Renderer::offsetP(const Eigen::Vector3d& p, const Eigen::Vector3d& n) const {
    const double eps = 1e-4; // シーンスケールに応じて
    return p + eps * n;
}

// ========== ここから “1点入射のラジアル・プロファイル検証” 本体 ==========
//
// ・半無限平板の 1 点 (xi,ni) から“正面入射”を仮定（=入射cos=1）。
// ・媒質内を等方散乱・NEEなし・Beer減衰で純ランダムウォーク。
// ・外に出たら、出射点 xo を記録。r=|xo-xi| のヒストグラムを作り Rd_emp を算出。
// ・理論値は evaluateBSSRDF(xi, xi+r*tangent) を呼んで Rd_theory として出す。
//
void Renderer::ProfileRadialBSSRDF_RW(const Eigen::Vector3d& xi,
                                      const Eigen::Vector3d& ni,
                                      const Material& mat,
                                      int bodyId,
                                      size_t samples,
                                      double rMax,
                                      int bins,
                                      const std::string& csvPath) const
{
    // 物性（スケール反映）
    const double scale = mat.scale;
    const double g=0.8;
    const Eigen::Vector3d sigma_a = mat.sigmas[0] * scale;   // RGB
    const Eigen::Vector3d sigma_s = (1-g)*mat.sigmas[1] * scale;   // RGB
    const Eigen::Vector3d sigma_t_v = sigma_a + sigma_s;     // RGB
    const double eta = mat.eta;

    // 代表値（サンプリング用）。RGBを厳密にやるなら追補可。
    const double sigma_t = std::max(1e-12, sigma_t_v.mean());
    const Eigen::Vector3d albedo_v = sigma_s.cwiseQuotient(sigma_t_v); // RGB
    const double albedo = std::clamp(albedo_v.mean(), 0.0, 0.999);

    const double dr = rMax / double(bins);
    std::vector<double> counts(bins, 0.0); // 重み付きカウント
    size_t exits = 0;

    // 入射側 Fresnel 透過（正面入射を仮定：cos=1）
    const double Ft_i = 1.0 - fresnelDielectric_R(1.0, eta, 1.0);

    for (size_t s = 0; s < samples; ++s) {
        // 表面直下へ入射（法線方向に入る）
        Eigen::Vector3d p = offsetP(xi, -ni);
        Eigen::Vector3d w = -ni; // 内部へ

        // スループット（RGB。初期は Ft_i）
        Eigen::Vector3d T(Ft_i, Ft_i, Ft_i);

        // ランダムウォーク
        for (int step = 0; step < 1<<20; ++step) {
            // 1) 自由行程
            const double t = sampleFreeFlight(urand(), sigma_t);
            Eigen::Vector3d p_next = p + t * w;

            // 2) 途中で境界に当たるかチェック
            {
                Ray r; r.org = p; r.dir = w;
                RayHit h;
                if (hitScene(r, h)) {
                    const double dSurf = (h.point - p).norm();
                    if (dSurf < t) {
                        // Beer（境界まで）
                        const Eigen::Vector3d Tr = ( (-sigma_a * dSurf).array().exp() ).matrix();
                        T = T.cwiseProduct(Tr);

                        // 内→外のフレネル
                        const Eigen::Vector3d n_b = h.normal.normalized();
                        const double cosI = std::max(0.0, -n_b.dot(w));
                        const double Rb   = fresnelDielectric_R(eta, 1.0, cosI);
                        const double Ft_o = 1.0 - Rb;

                        // 透過で外へ
                        if (urand() < Ft_o) {
                            Eigen::Vector3d wt;
                            if (refract(w, -n_b, eta, 1.0, wt)) {
                                // 出口点
                                const Eigen::Vector3d xo = h.point;
                                const double rlen = (xo - xi).norm();
                                if (rlen < rMax) {
                                    const int b = int(rlen / dr);
                                    if (b >= 0 && b < bins) {
                                        // 現在のRGBスループットを Luminance/平均で集計（どちらでもOK）
                                        const double wgt = (T[0] + T[1] + T[2]) / 3.0;
                                        counts[b] += wgt;
                                        exits++;
                                    }
                                }
                                // このフォトンは終了
                                break;
                            } else {
                                // 数値的に希。内部反射扱いへフォールバック
                            }
                        }

                        // 内部反射（または屈折失敗） → 反射方向へ
                        const Eigen::Vector3d w_ref = (w - 2.0 * w.dot(n_b) * n_b).normalized();
                        p = offsetP(h.point, -n_b);
                        w = w_ref;
                        // Beer は境界までで掛けた
                        continue;
                    }
                }
            }

            // 3) 境界前に媒質内で散乱
            const Eigen::Vector3d Tr = ( (-sigma_a * t).array().exp() ).matrix();
            T = T.cwiseProduct(Tr);

            // 散乱イベント：等方（Neumann展開の係数としてアルベドを掛ける近似）
            T = T.cwiseProduct(albedo_v);
            p = p_next;
            w = sampleIsotropicDir(urand(), urand());

            // ロシアンルーレット
            if (step > 8) {
                const double q = std::clamp(std::max({T[0],T[1],T[2]}), 0.05, 0.99);
                if (urand() > q) break;
                T /= q;
            }
        }
    }

    // --- 経験的 Rd(r) を作成（area normalization）---
    // Rd_emp(r_i) ≈ (Σ weights_in_bin / samples) / (2π r_i Δr)
    // 理論 Rd_theory(r_i) は evaluateBSSRDF で評価（RGB平均）
    std::ofstream ofs(csvPath);
    ofs << std::setprecision(10);
    ofs << "r_center,Rd_empirical,Rd_theory(Dipole),Rd_theory(Multipole),\n";

    // 接平面上の任意接線
    const Eigen::Vector3d tangent = ni.unitOrthogonal();

    for (int i = 0; i < bins; ++i) {
        const double r  = (i + 0.5) * dr;
        const double area = std::max(1e-16, (double)EIGEN_PI* 2* r * dr);
        const double Rd_emp = (counts[i] / double(samples)) / area;

        // 理論：xi から距離 r の点 xo_i（接平面上）で Rd を評価
        const Eigen::Vector3d xo_i = xi + r * tangent;
        const Color Rd_rgb = evaluateBSSRDF(xi, xo_i, mat);
        const Color Rd_rgb_m=evaluateBSSRDF_Multipole(xi, xo_i, mat,5,rMax,3);
        const double Rd_theory_Dipole = std::max(0.0, (Rd_rgb[0] + Rd_rgb[1] + Rd_rgb[2]) / 3.0);
        const double Rd_theory_Multipole = std::max(0.0, (Rd_rgb_m[0] + Rd_rgb_m[1] + Rd_rgb_m[2]) / 3.0);

        ofs << r << "," << Rd_emp << "," << Rd_theory_Dipole <<","<<Rd_theory_Multipole<< "\n";
    }
    ofs.close();

    std::cout << "[ProfileRadialBSSRDF_RW] samples=" << samples
              << " exits=" << exits
              << " csv=" << csvPath << std::endl;
}

#ifndef TWO_PI
#define TWO_PI (double)(2.0 * EIGEN_PI)
#endif




// 安全のためここで定義（ヘッダ不要）
static constexpr double TWO_PI_CNST = 6.28318530717958647692;


void Renderer::_ProfileRadialBSSRDF_RW(const Eigen::Vector3d& xi,
                                       const Eigen::Vector3d& ni,
                                       const Material& mat,
                                       int /*bodyId*/,
                                       size_t samples,
                                       double thickness,   // [m]
                                       double rMax,        // [m]
                                       int bins,
                                       const std::string& csvPath) const
{
    // --- 材質（単位を [1/m] に統一）---
    // 例：mat.sigmas が [1/mm] なら、mat.scale=1000 を入れておく。既に [1/m] なら 1.0。
    const double scale = mat.scale*1000; // ★ ここで *1000 はしない（mat 側で用意）
    const Eigen::Vector3d sigma_a_v = mat.sigmas[0] * scale;  // [1/m]
    const Eigen::Vector3d sigma_s_v = mat.sigmas[1] * scale;  // [1/m]（等方なら σs' = σs）
    const Eigen::Vector3d sigma_t_v = sigma_a_v + sigma_s_v;  // [1/m]
    const double sigma_t = std::max(1e-12, sigma_t_v.mean()); // サンプリング用代表値
    const double lt = 1.0 / sigma_t;                          // 平均自由行程 [m]
    const double eta = mat.eta;

    //std::cout << std::scientific
    //  << "[diag] sigma_a=" << sigma_a_v.transpose()
    //  << " sigma_s=" << sigma_s_v.transpose()
    //  << " sigma_t_mean=" << sigma_t
    //  << "  l_t=" << lt << " [m]"
    //  << "  T=" << thickness << " [m] (T/l_t=" << thickness/lt << ")"
    //  << std::endl;

    // --- ビン ---
    const double dr = rMax / std::max(1, bins);
    std::vector<double> counts(bins, 0.0);
    size_t exits = 0;

    // --- 基底（入射面の法線基準の局所座標） ---
    const Eigen::Vector3d n = ni.normalized();      // 外向き
    const Eigen::Vector3d u = n.unitOrthogonal();   // 接線1（tangent）
    auto depth = [&](const Eigen::Vector3d& P){ return (P - xi).dot(-n); }; // 入射面からの内向き距離 d>=0
    auto w_in  = [&](const Eigen::Vector3d& W){ return W.dot(-n); };        // 進行方向の内向き成分
    auto radial = [&](const Eigen::Vector3d& P){
        Eigen::Vector3d dP = P - xi;
    //debug dp=0になっているのか？
        //std::cout<<dP<<std::endl;

        Eigen::Vector3d tP = dP - n * dP.dot(n); // 平面へ正射影
        return tP.norm();
    };

    // 入射時フレネル透過（正面入射近似で十分）
    const double Ft_i = 1.0 - fresnelDielectric_R(1.0, eta, 1.0);

    // --- デバッグカウンタ（厚みが効いているかの確認用） ---
    size_t hitTop=0, hitBot=0, reflTop=0, reflBot=0, transTop=0, transBot=0;

    // --- 乱数トレーサ ---
    for (size_t s = 0; s < samples; ++s) {
        // 表面直下・内向きで開始
        Eigen::Vector3d p = xi + (-n) * 1e-4;
        Eigen::Vector3d w = -n;

        // RGBスループット（ここは平均で集計）
        Eigen::Vector3d T(Ft_i, Ft_i, Ft_i);

        for (;;) {
            // 自由行程は σt（reduced 媒質の作法：散乱は必ず起きる、吸収は区間 Beer のみ）
            double t = sampleFreeFlight(urand(), sigma_t);

            // 自由行程 t を境界でクリップしながら進める
            for (;;) {
                const double d  = depth(p);
                const double s_in = w_in(w);

                // 次に当たる境界まで（レイパラメータ）
                double dTop = std::numeric_limits<double>::infinity();
                double dBot = std::numeric_limits<double>::infinity();
                if (s_in < -1e-12) dTop =  d / (-s_in);              // 上面へ
                if (s_in >  1e-12) dBot = (thickness - d) / ( s_in); // 下面へ

                const bool hitTopBoundary = (dTop < dBot);
                const double dB = std::min(dTop, dBot);

                if (dB < t) {
                    // --- 先に境界に命中：Beer は境界まで ---
                    const Eigen::Vector3d Tr = ((-sigma_a_v * dB).array().exp()).matrix();
                    T = T.cwiseProduct(Tr);

                    // 境界まで進め、面上に厳密クランプ（接線成分は保持）
                    p += dB * w;
                    {
                        Eigen::Vector3d dP = p - xi;
                        Eigen::Vector3d tP = dP - n * dP.dot(n);
                        const double newDepth = hitTopBoundary ? 0.0 : thickness;
                        p = xi + tP + (-n) * newDepth;
                    }

                    // 境界の外向き法線
                    const Eigen::Vector3d nb = hitTopBoundary ? n : (-n);
                    if (hitTopBoundary) ++hitTop; else ++hitBot;

                    // Fresnel（内→外）。cos は nb に対する入射余弦そのまま（abs 不要）
                    const double cos_i = std::max(0.0, w.dot(nb));
                    const double Rb    = fresnelDielectric_R(eta, 1.0, cos_i);
                    const double Ft_o  = 1.0 - Rb;

                    if (urand() < Ft_o) {
                        // 透過 → 出射記録（フラックスなので cos を掛ける）
                        if (hitTopBoundary) ++transTop; else ++transBot;

                        const double cosNo = cos_i; // 屈折方向を明示しない簡易扱い
                        const double r = radial(p);

                        //debug countsが何故機能していないのか？

                        if (r < rMax) {
                            //std::cout<<"r="<<r<<":rMax="<<rMax<<std::endl;

                            const int b = int(r / dr);
                            if (b >= 0 && b < bins) {
                                //std::cout<<"r="<<r<<":rMax="<<rMax<<std::endl;


                                const double wgt = ((T[0]+T[1]+T[2]) * (1.0/3.0)) * cosNo;
                                counts[b] += wgt;
                            }
                        }
                        ++exits;
                        goto NEXT_PHOTON; // このフォトン終了
                    } else {
                        // 反射（鏡面反射：nb を使う）→ 残り距離で継続
                        if (hitTopBoundary) ++reflTop; else ++reflBot;
                        w = (w - 2.0 * w.dot(nb) * nb).normalized();
                        t -= dB;
                        if (t <= 1e-12) break; // 残距離ほぼ無し → 散乱へ
                        continue;               // 残距離で再度境界チェック
                    }
                } else {
                    // --- 境界に届かず媒質内で散乱 ---
                    const Eigen::Vector3d Tr = ((-sigma_a_v * t).array().exp()).matrix();
                    T = T.cwiseProduct(Tr);
                    p += t * w;

                    // reduced 媒質：離散の吸収ブランチは不要（Beer のみで吸収）
                    // 散乱は必ず発生
                    w = sampleIsotropicDir(urand(), urand());
                    break; // 次の自由行程へ
                }
            }

            // 検証ではロシアンルーレットは原則オフで OK
        }
        NEXT_PHOTON: ;
    }

    // デバッグ出力
    std::cout << "T="<<thickness
              << " hitTop="<<hitTop<<" hitBot="<<hitBot
              << " reflTop="<<reflTop<<" reflBot="<<reflBot
              << " transTop="<<transTop<<" transBot="<<transBot
              << std::endl;

    // --- CSV 出力：Rd_empirical と Multipole 理論の比較 ---
    std::ofstream ofs(csvPath);
    ofs << std::setprecision(10);
    ofs << "r_center,Rd_empirical,Rd_theory\n";



    const Eigen::Vector3d tangent = u; // 任意接線

    const double _dr=dr*1000;

    for (int i = 0; i < bins; ++i) {
        const double r = (i + 0.5) * _dr;
        const double area = std::max(1e-16, TWO_PI_CNST * r * _dr);
        const double Rd_emp = (counts[i] / double(samples)) / area;

        const double Fdr = Fdr_from_eta(mat.eta);
        const double C_boundary = 1-Fdr ; // 入出 2 回


        // 理論（Multipole）。★ thickness は [m] のまま渡す（*1000 しない）
        const Eigen::Vector3d xo = xi + r * tangent;
        const Color Rd_th_c = evaluateBSSRDF_Multipole(xi, xo, mat, thickness*1000, rMax,/*M=*/3);
        const double Rd_theory = (Rd_th_c[0] + Rd_th_c[1] + Rd_th_c[2]) *C_boundary* (1.0/3.0);


        ofs << r << "," << Rd_emp << "," << Rd_theory << "\n";
    }
    ofs.close();

    std::cout << "[_ProfileRadialBSSRDF_RW] samples=" << samples
              << " exits=" << exits
              << " thickness=" << thickness
              << " csv=" << csvPath << std::endl;
}
// Rd(r) の放射プロファイル（方針A: Fresnelは外で別に掛ける）
double Renderer::computeRhoBSSRDF(const Eigen::Vector3d& xi, const Material& mat,double rMax, int N, bool includeFresnelOutside) const{
    const double dr = rMax / N;
    auto Rd_scalar = [&](double r)->double {
        double Rd = evaluateBSSRDFScalarAtDistance(r, mat);                     // RGB → スカラー（平均など）
        return Rd;
    };

    double integral = 0.0;
    double prev = Rd_scalar(0.0);
    for (int i = 1; i <= N; ++i) {
        double r = i * dr;
        double cur = Rd_scalar(r);
        double a = 2.0 * EIGEN_PI * (r - dr) * prev;
        double b = 2.0 * EIGEN_PI * r * cur;
        integral += 0.5 * (a + b) * dr;  // 台形公式で ∫ 2π r Rd(r) dr
        prev = cur;
    }

    if (includeFresnelOutside) {
        double Ft_in  = Fdr_from_eta( mat.eta);  // 半球平均の近似

        integral *= (Ft_in * Ft_in);
    }
    return integral; // これが ρ_BSSRDF
}



void  Renderer::Kensyou(const unsigned int &samples) const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());
    int p_x=image.width/2;
    int p_y=image.height/2;
    Ray ray;
    RayHit hit;
    camera.filmView(p_x, p_y, ray);

    if (hitScene(ray, hit)) {
        const Eigen::Vector3d xi = hit.point;
 //
 //      ProfileRadialBSSRDF_RW(xi,hit.normal,bodies[hit.idx].material,hit.idx,
 //       /*samples=*/samples,  // 100万くらいから（数分〜）
//
//    /*rMax=*/20,          // 2cm相当（材質/スケールに合わせて）
//    /*bins=*/300,
//    /*csvPath=*/"kensyou_4.csv");

        const double dr = 0.02 / std::max(1, 300);
        std::ofstream ofs("kensyou_7.csv");
        ofs << std::setprecision(10);
        ofs << "r_center,Dipole,Rd_theory_M\n";
        const Eigen::Vector3d n = hit.normal.normalized();      // 外向き
        const Eigen::Vector3d u = n.unitOrthogonal();   // 接線1（tangent）
        const Eigen::Vector3d tangent = u; // 任意接線

        const double _dr=dr*1000;

        for (int i = 0; i < 300; ++i) {
            const double r = (i + 0.5) * _dr;
            const double area = std::max(1e-16, TWO_PI_CNST * r * _dr);


            const double Fdr = Fdr_from_eta(bodies[hit.idx].material.eta);
            const double C_boundary = 1-Fdr ; // 入出 2 回


            // 理論（Multipole）。★ thickness は [m] のまま渡す（*1000 しない）
            const Eigen::Vector3d xo = xi + r * tangent;
            const double Rd_theory =evaluateBSSRDFScalarAtDistance(r,bodies[hit.idx].material);
            const Color Rd_th_c = evaluateBSSRDF_Multipole(xi, xo, bodies[hit.idx].material, 5*1000, 0.02,/*M=*/3);
            const double Rd_theory_M = (Rd_th_c[0] + Rd_th_c[1] + Rd_th_c[2]) *C_boundary* (1.0/3.0);


            ofs << r << "," << Rd_theory<< "," << Rd_theory_M << "\n";
        }
        ofs.close();
        double ppp=computeRhoBSSRDF(xi,bodies[hit.idx].material,0.02,samples*10,false);
        Color kddddd=bodies[hit.idx].getKd();
        double qqq=(kddddd[0]+kddddd[1]+kddddd[2])/3;
        std::cout<<"kd from mat="<<qqq<<", p_bssrdf="<<ppp<<std::endl;

        double thickness = 0.005;  // 2mm スラブ
        _ProfileRadialBSSRDF_RW(
            xi, hit.normal,bodies[hit.idx].material,hit.idx,
            /*samples=*/samples*10,
            /*thickness=*/thickness,
            /*rMax=*/0.02,
            /*bins=*/300,
            /*csvPath=*/"kensyou_5.csv"
        );

        _ProfileRadialBSSRDF_RW(
            xi, hit.normal,bodies[hit.idx].material,hit.idx,
            /*samples=*/samples*10,
            /*thickness=*/0.002,
            /*rMax=*/0.02,
            /*bins=*/300,
            /*csvPath=*/"kensyou_6_1.csv"
        );
        _ProfileRadialBSSRDF_RW(
            xi, hit.normal,bodies[hit.idx].material,hit.idx,
            /*samples=*/samples*10,
            /*thickness=*/0.005,
            /*rMax=*/0.02,
            /*bins=*/300,
            /*csvPath=*/"kensyou_6_2.csv"
        );
        _ProfileRadialBSSRDF_RW(
            xi, hit.normal,bodies[hit.idx].material,hit.idx,
            /*samples=*/samples*10,
            /*thickness=*/0.02,
            /*rMax=*/0.02,
            /*bins=*/300,
            /*csvPath=*/"kensyou_6_3.csv"
        );






    }



}


void Renderer::computeLocalFrame(const Eigen::Vector3d &w, Eigen::Vector3d &u, Eigen::Vector3d &v) {
    if(fabs(w.x()) > 1e-3)
        u = Eigen::Vector3d::UnitY().cross(w).normalized();
    else
        u = Eigen::Vector3d::UnitX().cross(w).normalized();

    v = w.cross(u);
}
