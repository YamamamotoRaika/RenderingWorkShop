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

    loadTex("../Day3/lpshead/face_beard_density_test.jpg");
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
bool Renderer::_hitScene(const Ray &ray, RayHit &hit) const {
    /// hitするBodyのうち最小距離のものを探す
    hit.t = DBL_MAX;
    hit.idx = -1;
    for(int i = 0; i < bodies.size();++ i) {
        if (bodies[i].isBeard()) continue;
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
#pragma omp parallel for
    for (int p_y = 0; p_y < image.height; ++p_y) {
        for (int p_x = 0; p_x < image.width; ++p_x) {
            const int p_idx = p_y * image.width + p_x;
            Ray ray; RayHit hit;
            camera.filmView(p_x, p_y, ray);

            if (hitScene(ray, hit)) {
                Color reflectRadiance = Color::Zero();
                const Color kd = bodies[hit.idx].getKd();
                for (int i = 0; i < samples; ++i) {
                    Ray _ray; RayHit _hit;
                    diffuseSample(hit.point, hit.normal, _ray);
                    if (hitScene(_ray, _hit)) {
                        reflectRadiance += kd.cwiseProduct(bodies[_hit.idx].getEmission());
                    }
                }
                image.pixels[p_idx] = bodies[hit.idx].getEmission()
                    + reflectRadiance / static_cast<double>(samples);
            } else {
                image.pixels[p_idx] = bgColor;
            }
        }
    }
    return image;
}
inline Eigen::Vector3d d(Eigen::Vector3d x,Eigen::Vector3d org ) {
    Eigen::Vector3d d=x-org;;


    return d;
}

Image Renderer::_directIlluminationRender(const unsigned int &samples) const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());



    /// フィルム上のピクセル全てに向けてレイを飛ばす
#pragma omp parallel for
    for(int p_y = 0; p_y < image.height; p_y++) {
        int Count_x=0;
        int Count_xx=0;
        int Count_xxx=0;

        for(int p_x = 0; p_x < image.width; p_x++) {
            const int p_idx = p_y * image.width + p_x;
            Ray ray;
            RayHit hit;
            camera.filmView(p_x, p_y, ray);

            if (hitScene(ray, hit)) {
                //debug you



                if(bodies[hit.idx].isLight()) {
                    image.pixels[p_idx] = bodies[hit.idx].getEmission();
                } else {
                    Color reflectRadiance = Color::Zero();
                    for (int i = 0; i < samples; ++i) {
                        /// 衝突点hit.pointから半球上のランダムな方向にレイを飛ばす
                        Ray _ray; RayHit _hit;
                        constexpr double eps = 1e-4; // シーンスケールに合わせて 1e-4〜1e-3 で調整
                        Eigen::Vector3d org = hit.point + hit.normal ;
                        //std::cout<<hit.normal<<std::endl;




                        diffuseSample(hit.point, hit.normal, _ray); // hit.point じゃなく org を使う


                        /// もしBodyに当たったら,その発光量を加算する
                        if (hitScene(_ray, _hit) ) {
                            //kokowomiru

                            Count_x++;
                            if (hit.idx==_hit.idx) {

                                if (bodies[hit.idx].getId()!=0&&bodies[hit.idx].getId()==bodies[_hit.idx].getId()) {
                                    Count_xxx++;
                                }
                            }
                            if (bodies[_hit.idx].isLight())  Count_xx++;

                            //std::cout<< bodies[hit.idx].getKd()<<std::endl;
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
        std::cout << "***"<<Count_x << " ," << Count_xx << " ," << Count_xxx << "***"<<std::endl;
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

void Renderer::SSSSample(const Eigen::Vector3d &incidentPoint,
                         const Eigen::Vector3d &normal,
                         const Material& material,
                         RayHit hit,
                         Ray &out_Ray,double &out_pdfA) const {
    const auto& tbl = getRadialTable(material);

    // ---- 1. 半径と方位のサンプル（出口点の位置） ----
    const double uR   = (double)rand() / (double)RAND_MAX;
    const double uPhi = rand() ;

    const double r   = tbl.sampleR(uR);              // 半径
    const double phi = 2.0 * EIGEN_PI * uPhi;        // 平面上の角度

    Eigen::Vector3d u, v;
    computeLocalFrame(normal, u, v);                 // 入射法線まわりの接平面

    const double dx = r * std::cos(phi);
    const double dy = r * std::sin(phi);

    // 接平面上の候補点
    Eigen::Vector3d pCandidate = incidentPoint + dx * u + dy * v;

    // ---- 2. 同じBodyの実表面に投影 ----
    Eigen::Vector3d xo, n_o;
    if (!projectToSurface(pCandidate, normal, hit.idx, xo, n_o)) {
        // 投影に失敗したらこのサンプルは寄与ゼロ扱い
        out_Ray.org  = incidentPoint;
        out_Ray.dir  = Eigen::Vector3d::Zero();
        out_pdfA     = 0.0;
        return;
    }

    // ---- 3. 出口点での方向サンプル（コサイン半球） ----
    const double uDir1 = (double)rand() / (double)RAND_MAX;
    const double uDir2 = rand() ;

    const double phiDir   = 2.0 * EIGEN_PI * uDir1;
    const double cosTheta = std::sqrt(1.0 - uDir2);
    const double sinTheta = std::sqrt(uDir2);

    Eigen::Vector3d du, dv;
    computeLocalFrame(n_o, du, dv);                  // ★ 出口法線 n_o で基底を作る

    const Eigen::Vector3d wi =
        sinTheta * std::cos(phiDir) * du +
        cosTheta * n_o +
        sinTheta * std::sin(phiDir) * dv;

    const double eps = 1e-4;
    out_Ray.org = xo + n_o * eps;                    // 自己交差防止
    out_Ray.dir = wi.normalized();

    // ---- 4. 面積PDF p_A(r) の計算 ----
    // Multipole BSSRDF を使って Rd(r) を評価
    /*
    const Color RdC = evaluateBSSRDF_Multipole(incidentPoint, xo, material,
                                               200,
                                               0.2,
                                              5.0);
    */
    const Color RdC =evaluateBSSRDF(incidentPoint,out_Ray.org,material);
    const double Rd = (RdC[0] + RdC[1] + RdC[2]) / 3.0;
    //std::cout <<"Rd=" <<Rd << std::endl;

    // tbl.Z を ∫_0^∞ 2π r Rd(r) dr として定義している前提なら:
    const double pA = 2.0 * EIGEN_PI * r * Rd / tbl.Z;
    out_pdfA = pA;

    // 方向PDF（cosine）は cosTheta / π なので、
    // それは呼び出し側で別途使う or ここで一緒に返す設計でもOK

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
    //projectToSurface(,,,out_Ray.org,)
    // ★ 面積pdf： p_A = Rd(r) / Z
    const double Rd = evaluateBSSRDFScalarAtDistance(r,material);
   // const Color RdC = evaluateBSSRDF_Multipole(incidentPoint,out_Ray.org,material,200,0.2,5);
   // const double Rd=(RdC[1]+RdC[2]+RdC[0])/3.0;

    out_pdfA = std::max(1e-18, Rd / tbl.Z);

    //out_pdfA = std::max(1e-18, Rd / tbl.Z);
}

void Renderer::_Multipole_kensyou(const Eigen::Vector3d &incidentPoint,
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

    std::vector<SSSLayer> skin = {
        // 表皮（薄い・吸収強め）
        { {0.44,0.22,0.14}, {1.2,1.3,1.4}, 0.8, 1.4, 0.2 },
        // 真皮（厚い・散乱強）
        { {0.02,0.03,0.04}, {2.5,3.0,3.5}, 0.8, 1.4, 3.0 }
    };

    MLProfileSettings opt;
    opt.imagesPerSide = 8;

    double rr = 0.5; // mm
    Eigen::Vector3d RdC = Rd_multipole_multilayer(skin, rr, opt);



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
    /*
    std::cout<<"eta="<<eta<<std::endl;
    std::cout<<"scale="<<scale<<std::endl;
    std::cout<<"sigma_s="<<sigma_s<<std::endl;
    std::cout<<"sigma_a="<<sigma_a<<std::endl;
*/
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
    //std::cout << "Rd:"<< Rd[0]<<":"<<Rd[1]<<":"<<Rd[2]<< std::endl;
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


//なんか壊れてるわこいつ
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
                        const Color kd = bodies[hit.idx].getKd();

                        for (int i = 0; i < samples; ++i) {
                        // … ピクセル内ループ …
                        Ray _ray; RayHit _hit;

                        // ① 位置：Rd に沿って
                        double pdfA = 0.0;
                        _SSSSample(hit.point, hit.normal, bodies[hit.idx].getMaterial(), _ray, pdfA); // xi = _ray.org
                        //SSSSample(hit.point, hit.normal, bodies[hit.idx].getMaterial(),hit, _ray, pdfA); // xi = _ray.org

                        // ② 方向：既存の diffuseSample を使う（コサイン分布）
                        //    → 後で p_dir = cosθ/π で割る
                        diffuseSample(_ray.org, hit.normal, _ray);
                        const double cos_i = std::max(0.0, hit.normal.dot(_ray.dir));



                        if (hitScene(_ray, _hit) && bodies[_hit.idx].isLight()) {
                            const Color Li = bodies[_hit.idx].getEmission();

                            // ③ S = (Ft_i Ft_o / π) * Rd_rgb(d)
                            const double cos_o = std::max(0.0, hit.normal.dot(-ray.dir)); // カメラ側
                            const double Ft_i  = fresnelT_schlick(bodies[hit.idx].getMaterial().eta, cos_i);
                            const double Ft_o  = fresnelT_schlick(bodies[hit.idx].getMaterial().eta, cos_o);


                            const Color Rd_rgb = evaluateBSSRDF(_ray.org, hit.point, bodies[hit.idx].getMaterial());
                            const Color S_full = (Ft_i * Ft_o / EIGEN_PI) * Rd_rgb;

                            //std::cout << (Ft_i * Ft_o / EIGEN_PI)<<std::endl;


                            // ④ 不偏推定：
                            reflectRadiance += (Ft_i * Ft_o ) * kd.cwiseProduct(Li) / (pdfA);
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
                        const Color kd = bodies[hit.idx].getKd();

                        for (int i = 0; i < samples; ++i) {
                        Ray _ray; RayHit _hit;

                        // ★ 基準点は hit.point を渡す（未初期化の _ray.org ではない）
                        double pdfA = 0.0;

                        // ③ S = (Ft_i Ft_o / π) * Rd_rgb(d)
                        const double cos_i = std::max(0.0, hit.normal.dot(_ray.dir));
                        const double cos_o = std::max(0.0, hit.normal.dot(-ray.dir)); // カメラ側
                        const double Ft_i  = fresnelT_schlick(bodies[hit.idx].getMaterial().eta, cos_i);
                        const double Ft_o  = fresnelT_schlick(bodies[hit.idx].getMaterial().eta, cos_o);

                        //マルチポール複数層検証用
                        //_Multipole_kensyou(hit.point, hit.normal, bodies[hit.idx].material, _ray, pdfA);
                        //SSSSample(hit.point, hit.normal, bodies[hit.idx].getMaterial(),hit, _ray, pdfA);
                        _SSSSample(hit.point, hit.normal, bodies[hit.idx].getMaterial(), _ray, pdfA); // xi = _ray.org

                        // 出射方向のサンプル（従来通り）
                        diffuseSample(_ray.org, hit.normal, _ray);



                        if (hitScene(_ray, _hit) && bodies[_hit.idx].isLight()) {
                            const Color Li = bodies[_hit.idx].getEmission();                 // 光源放射
                            const Color S  = evaluateBSSRDF(_ray.org, hit.point, bodies[hit.idx].material); // Rd評価（RGB）
                            // ★ 面積pdfで割る
                           // reflectRadiance += bodies[hit.idx].getKd().cwiseProduct(S.cwiseProduct(Li)) / pdfA;
                            reflectRadiance += (Ft_i*Ft_o)*bodies[hit.idx].getKd().cwiseProduct(Li) / pdfA ;
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

Image Renderer::KAI_SSSdirectIlluminationRender(const unsigned int &samples) const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());
    //loadTex("beardTex.jpg");

    Material  beardMat=Material(
    Color(0.1, 0.1, 0.4),
    0.8,                            // kd（とりあえず肌と同じ）
    0.0,                            // emission なし
    1.55,                           // eta: 髪/毛の屈折率っぽく少し高め
    // sigma_a: 吸収係数（肌よりかなり大きく、特にB成分を強く吸う）
    Eigen::Vector3d(0.40, 0.55, 0.40),
    // sigma_s': 拡散係数（肌より小さめで、あまり遠くまで拡散しない）
    Eigen::Vector3d(0.05, 0.10, 0.30),
    true                            // isSubsurface
);

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
                        const Color kd = bodies[hit.idx].getKd();

                        for (int i = 0; i < samples; ++i) {
                        // … ピクセル内ループ …
                        Ray _ray; RayHit _hit;

                        // ① 位置：Rd に沿って
                        double pdfA = 0.0;
                        bool isBeard=true;
                       // _SSSSample(hit.point, hit.normal, bodies[hit.idx].getMaterial(), _ray, pdfA); // xi = _ray.org
                       //double beardWeight =beardDensity(bodies[hit.idx].getTriangle().getUV(hit.point));
                        double beardWeight =beardDensity(bodies[hit.idx].getTriangle().getUV(hit.point));
                       //double beardWeight =0;

                        //std::cout<<"hit_point="<<hit.point<<std::endl;
                        //std::cout<<"beardWeight="<<beardWeight<<std::endl;

                       SSSSampleSkinBeardMixed(hit.point, hit.normal, bodies[hit.idx].getMaterial(),beardMat,beardWeight, _ray, pdfA,isBeard);
                        // ② 方向：既存の diffuseSample を使う（コサイン分布）
                        //    → 後で p_dir = cosθ/π で割る
                        diffuseSample(_ray.org, hit.normal, _ray);
                        const double cos_i = std::max(0.0, hit.normal.dot(_ray.dir));



                        if (_hitScene(_ray, _hit) && bodies[_hit.idx].isLight()) {
                            const Color Li = bodies[_hit.idx].getEmission();

                            // ③ S = (Ft_i Ft_o / π) * Rd_rgb(d)
                            const double cos_o = std::max(0.0, hit.normal.dot(-ray.dir)); // カメラ側
                            const double Ft_i  = fresnelT_schlick(bodies[hit.idx].getMaterial().eta, cos_i);
                            const double Ft_o  = fresnelT_schlick(bodies[hit.idx].getMaterial().eta, cos_o);


                            const Color Rd_rgb = evaluateBSSRDF(_ray.org, hit.point, bodies[hit.idx].getMaterial());
                            const Color S_full = (Ft_i * Ft_o / EIGEN_PI) * Rd_rgb;

                            //std::cout << (Ft_i * Ft_o / EIGEN_PI)<<std::endl;


                            // ④ 不偏推定
                            reflectRadiance += (Ft_i * Ft_o )*kd.cwiseProduct( Li ) / (pdfA );
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


Image Renderer::_KAI_SSSdirectIlluminationRender(const unsigned int &samples) const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());
    //loadTex("beardTex.jpg");

    Material  beardMat=Material(
    Color(0.1, 0.1, 0.4),
    0.8,                            // kd（とりあえず肌と同じ）
    0.0,                            // emission なし
    1.55,                           // eta: 髪/毛の屈折率っぽく少し高め
    // sigma_a: 吸収係数（肌よりかなり大きく、特にB成分を強く吸う）
    Eigen::Vector3d(0.40, 0.55, 0.40),
    // sigma_s': 拡散係数（肌より小さめで、あまり遠くまで拡散しない）
    Eigen::Vector3d(0.05, 0.10, 0.30),
    true                            // isSubsurface
);

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
                        Color kd = bodies[hit.idx].getKd();
                        if (bodies[hit.idx].type == ShapeType::Mesh) {
                            const Triangle tri = bodies[hit.idx].getTriangle();
                            const Eigen::Vector2d uv = tri.getUV(hit.point);
                            kd = bodies[hit.idx].mesh.getTexture(uv);
                        }


                        for (int i = 0; i < samples; ++i) {
                        // … ピクセル内ループ …
                        Ray _ray; RayHit _hit;

                        // ① 位置：Rd に沿って
                        double pdfA = 0.0;
                        bool isBeard=true;
                       // _SSSSample(hit.point, hit.normal, bodies[hit.idx].getMaterial(), _ray, pdfA); // xi = _ray.org
                       //double beardWeight =beardDensity(bodies[hit.idx].getTriangle().getUV(hit.point));
                        double beardWeight =beardDensity(bodies[hit.idx].getTriangle().getUV(hit.point));
                       //double beardWeight =0;

                        //std::cout<<"hit_point="<<hit.point<<std::endl;
                        //std::cout<<"beardWeight="<<beardWeight<<std::endl;

                       SSSSampleSkinBeardMixed(hit.point, hit.normal, bodies[hit.idx].getMaterial(),beardMat,beardWeight, _ray, pdfA,isBeard);
                        // ② 方向：既存の diffuseSample を使う（コサイン分布）
                        //    → 後で p_dir = cosθ/π で割る
                        diffuseSample(_ray.org, hit.normal, _ray);
                        const double cos_i = std::max(0.0, hit.normal.dot(_ray.dir));



                        if (_hitScene(_ray, _hit) && bodies[_hit.idx].isLight()) {
                            const Color Li = bodies[_hit.idx].getEmission();

                            // ③ S = (Ft_i Ft_o / π) * Rd_rgb(d)
                            const double cos_o = std::max(0.0, hit.normal.dot(-ray.dir)); // カメラ側
                            const double Ft_i  = fresnelT_schlick(bodies[hit.idx].getMaterial().eta, cos_i);
                            const double Ft_o  = fresnelT_schlick(bodies[hit.idx].getMaterial().eta, cos_o);


                            const Color Rd_rgb = evaluateBSSRDF(_ray.org, hit.point, bodies[hit.idx].getMaterial());
                            const Color S_full = (Ft_i * Ft_o / EIGEN_PI) * Rd_rgb;

                            //std::cout << (Ft_i * Ft_o / EIGEN_PI)<<std::endl;


                            // ④ 不偏推定
                            reflectRadiance += (Ft_i * Ft_o )*kd.cwiseProduct( Li ) / (pdfA );
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
bool Renderer::ProfileRadialBSSRDF_BeardMixed_CameraCenter(size_t angularSamples,
                                                           double rMax,
                                                           int bins,
                                                           const std::string& csvPath) const {
    if (angularSamples == 0 || bins <= 0 || rMax <= 0.0) {
        std::cout << "[RdProfile] invalid params\n";
        return false;
    }

    const int p_x = camera.getFilm().resolution.x() / 2;
    const int p_y = camera.getFilm().resolution.y() / 2;

    Ray ray;
    RayHit hit;
    camera.filmView(p_x, p_y, ray);
    if (!hitScene(ray, hit)) {
        std::cout << "[RdProfile] no hit at camera center\n";
        return false;
    }
    if (bodies[hit.idx].isLight()) {
        std::cout << "[RdProfile] camera center hit light\n";
        return false;
    }

    const Material skinMat = bodies[hit.idx].getMaterial();
    Material beardMat = Material(
        Color(0.1, 0.1, 1.6),
        0.8,
        0.0,
        1.55,
        Eigen::Vector3d(0.35, 0.55, 1.6),
        Eigen::Vector3d(0.05, 0.10, 0.15),
        true
    );

    if (bodies[hit.idx].type != ShapeType::Mesh) {
        std::cout << "[RdProfile] target body is not mesh\n";
        return false;
    }

    const Mesh& mesh = bodies[hit.idx].mesh;
    if (mesh.triangles.empty()) {
        std::cout << "[RdProfile] mesh has no triangles\n";
        return false;
    }

    const double Do = beardDensity(bodies[hit.idx].getTriangle().getUV(hit.point));

    std::ofstream ofs(csvPath);
    if (!ofs) {
        std::cout << "[RdProfile] failed to open " << csvPath << "\n";
        return false;
    }
    ofs << "r,Rd_skin,Rd_beard,Rd_mix\n";

    auto avgRGB = [](const Color& c) {
        return (c[0] + c[1] + c[2]) / 3.0;
    };

    std::vector<double> sum_skin(bins, 0.0);
    std::vector<double> sum_beard(bins, 0.0);
    std::vector<double> sum_mix(bins, 0.0);
    std::vector<size_t> count(bins, 0);

    std::vector<double> cdf(mesh.triangles.size(), 0.0);
    double areaSum = 0.0;
    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        areaSum += mesh.triangles[i].area();
        cdf[i] = areaSum;
    }
    if (areaSum <= 0.0) {
        std::cout << "[RdProfile] mesh area is zero\n";
        return false;
    }

    auto sampleTriangle = [&](double u)->const Triangle& {
        double target = u * areaSum;
        auto it = std::lower_bound(cdf.begin(), cdf.end(), target);
        size_t idx = std::min<size_t>(cdf.size() - 1, (size_t)std::distance(cdf.begin(), it));
        return mesh.triangles[idx];
    };

    auto samplePointOnTri = [&](const Triangle& tri, double u1, double u2)->Eigen::Vector3d {
        double su = std::sqrt(u1);
        double b0 = 1.0 - su;
        double b1 = su * (1.0 - u2);
        double b2 = su * u2;
        return tri.v0 * b0 + tri.v1 * b1 + tri.v2 * b2;
    };

    for (size_t s = 0; s < angularSamples; ++s) {
        const Triangle& tri = sampleTriangle(urand());
        const Eigen::Vector3d xo = samplePointOnTri(tri, urand(), urand());
        const double r = (xo - hit.point).norm();
        
        if (r <= 0.0 || r > rMax) continue;

        const double Di = beardDensity(tri.getUV(xo));
        double w = 1.0 - (1.0 - Do) * (1.0 - Di);
        w = std::clamp(w, 0.0, 1.0);

        const Color Rd_skin = evaluateBSSRDF(xo, hit.point, skinMat);
        const Color Rd_beard = evaluateBSSRDF(xo, hit.point, beardMat);
        const Color Rd_mix = (1.0 - w) * Rd_skin + w * Rd_beard;

        const int bin = std::min(bins - 1, std::max(0, int((r / rMax) * bins)));
        sum_skin[bin] += avgRGB(Rd_skin);
        sum_beard[bin] += avgRGB(Rd_beard);
        sum_mix[bin] += avgRGB(Rd_mix);
        count[bin] += 1;
         for (size_t s = 0; s < angularSamples; ++s) {
        const Triangle& tri = sampleTriangle(urand());
        const Eigen::Vector3d xo = samplePointOnTri(tri, urand(), urand());
        const double r = (xo - hit.point).norm();
        
        if (r <= 0.0 || r > rMax) continue;

        const double Di = beardDensity(tri.getUV(xo));
        double w = 1.0 - (1.0 - Do) * (1.0 - Di);
        w = std::clamp(w, 0.0, 1.0);

        const Color Rd_skin = evaluateBSSRDF(xo, hit.point, skinMat);
        const Color Rd_beard = evaluateBSSRDF(xo, hit.point, beardMat);
        const Color Rd_mix = (1.0 - w) * Rd_skin + w * Rd_beard;

        const int bin = std::min(bins - 1, std::max(0, int((r / rMax) * bins)));
        sum_skin[bin] += avgRGB(Rd_skin);
        sum_beard[bin] += avgRGB(Rd_beard);
        sum_mix[bin] += avgRGB(Rd_mix);
        count[bin] += 1;
         
    }

        
    }

    for (int i = 0; i < bins; ++i) {
        const double r = rMax * (double(i) + 0.5) / double(bins);
        if (count[i] == 0) {
            ofs << r << ",0,0,0\n";
           
        } else {
            ofs << r << ","
                << (sum_skin[i] / double(count[i])) << ","
                << (sum_beard[i] / double(count[i])) << ","
                << (sum_mix[i] / double(count[i])) << "\n";
        }
    }

    std::cout << "[RdProfile] saved " << csvPath << "\n";
    return true;
}

void Renderer::SSSSampleSkinBeardMixed(const Eigen::Vector3d& incidentPoint,
                                       const Eigen::Vector3d& normal,
                                       const Material& skinMat,
                                       const Material& beardMat,
                                       double beardWeight,
                                       Ray& out_Ray,
                                       double& out_pdfA,
                                       bool& out_isBeard) const
{
    // --- 0. カバレッジからミックス比を決定 ---
    // beardWeight: その点で髭が覆っている面積割合 (0 = 髭無し, 1 = 髭だけ)
    beardWeight = std::clamp(beardWeight, 0.0, 1.0);
    const double wBeard = beardWeight;
    const double wSkin  = 1.0 - beardWeight;

    // --- 1. 肌・髭それぞれのラジアルテーブルを取得 ---
    const auto& tblSkin  = getRadialTable(skinMat);
    const auto& tblBeard = getRadialTable(beardMat);

    const double Zs = tblSkin.Z;
    const double Zb = tblBeard.Z;

    // --- 2. どちらの成分からサンプルするか決める ---
    const double uComp = rand()/RAND_MAX;
    const bool chooseBeard = (uComp < wBeard);
    const BSSRDFRadialCDF& tbl = chooseBeard ? tblBeard : tblSkin;

    // --- 3. r, φ をサンプルして出口位置を決める ---
    const double u1  = (double)rand() / (double)RAND_MAX;
    const double u2  = (double)rand() / (double)RAND_MAX;
    const double r   = tbl.sampleR(u1);
    const double phi = 2.0 * EIGEN_PI * u2;

    Eigen::Vector3d t, b;
    computeLocalFrame(normal, t, b);

    const double dx = r * std::cos(phi);
    const double dy = r * std::sin(phi);

    out_Ray.org = incidentPoint + dx * t + dy * b;

    // --- 4. 面積 PDF を計算（混合分布の PDF を返す） ---
    //   p_A_skin(r)  = 2π r Rd_skin(r)  / Zs
    //   p_A_beard(r) = 2π r Rd_beard(r) / Zb
    //   p_mix(r) = wSkin * p_A_skin + wBeard * p_A_beard

    const double Rd_s = Renderer::evaluateBSSRDFScalarAtDistance(r, skinMat);
    const double Rd_b = Renderer::evaluateBSSRDFScalarAtDistance(r, beardMat);

    double pdfA_skin  = 0.0;
    double pdfA_beard = 0.0;
    if (Zs > 0.0) pdfA_skin  = std::max(1e-18, Rd_s / Zs);;
    if (Zb > 0.0) pdfA_beard = std::max(1e-18, Rd_b / Zb);;

    const double pdfMix = wSkin * pdfA_skin + wBeard * pdfA_beard;

    out_pdfA   = std::max(1e-18, pdfMix);
    out_isBeard = chooseBeard;
}
void Renderer::_SSSSampleSkinBeardMixed(
    const Eigen::Vector3d& incidentPoint,
    const Eigen::Vector3d& normal,
    const Material& skinMat,
    const Material& beardMat,
    double beardWeight,      // ρ_beard(x): 0〜1
    Ray& out_Ray,
    double& out_pdfA,
    bool& out_isBeard) const
{
    // 0. 髭密度をクランプ
    const double rho = std::clamp(beardWeight, 0.0, 1.0);
    const double wSkin  = 1.0 - rho;
    const double wBeard = rho;

    // 1. パラメータを線形補間して「混合マテリアル」を作る
    Material mixed = skinMat;   // ベースは肌マテリアルをコピー

    // 吸収・散乱係数を補間
    mixed.sigmas[0] = wSkin * skinMat.sigmas[0] + wBeard * beardMat.sigmas[0];
    mixed.sigmas[1] = wSkin * skinMat.sigmas[1] + wBeard * beardMat.sigmas[1];

    // 必要に応じて η や color も補間しても良い（好み）
    mixed.eta = wSkin * skinMat.eta + wBeard * beardMat.eta;
    mixed.color = wSkin * skinMat.color + wBeard * beardMat.color;

    // 2. 混合マテリアル用のラジアルテーブルを取得
    const BSSRDFRadialCDF& tbl = getRadialTable(mixed);
    const double Z = tbl.Z;

    // 3. r, φ をサンプリング
    const double u1 = (double)rand() / (double)RAND_MAX;
    const double u2 = (double)rand() / (double)RAND_MAX;
    const double r   = tbl.sampleR(u1);
    const double phi = 2.0 * EIGEN_PI * u2;

    Eigen::Vector3d t, b;
    computeLocalFrame(normal, t, b);

    const double dx = r * std::cos(phi);
    const double dy = r * std::sin(phi);

    out_Ray.org = incidentPoint + dx * t + dy * b;

    // 4. 面積 pdf を計算（単一成分なのでシンプル）
    const double Rd_mix = Renderer::evaluateBSSRDFScalarAtDistance(r, mixed);

    double pdfA = 0.0;
    if (Z > 0.0) {
        // 本来は 2πr を含めた形だけど、テーブルの定義に合わせてここは調整してね
        pdfA = std::max(1e-18, Rd_mix / Z);
    }

    out_pdfA = pdfA;
}


namespace {
    // [edge0, edge1] で 0→1 にスムーズに変化する smoothstep
    inline double smoothstep(double edge0, double edge1, double x)
    {
        double t = (x - edge0) / (edge1 - edge0);
        t = std::clamp(t, 0.0, 1.0);
        return t * t * (3.0 - 2.0 * t);
    }
} // anonymous namespace

double Renderer::beardDensity(const Eigen::Vector2d& uv) const
{
    if (m_beardW == 0 || m_beardH == 0) return 0.0;

    double u = uv.x() - std::floor(uv.x());
    double v = uv.y() - std::floor(uv.y());

    double x = u * (m_beardW  - 1);
    double y = (1.0 - v) * (m_beardH - 1);

    int x0 = std::clamp((int)std::floor(x), 0, m_beardW  - 1);
    int y0 = std::clamp((int)std::floor(y), 0, m_beardH - 1);
    int x1 = std::clamp(x0 + 1, 0, m_beardW  - 1);
    int y1 = std::clamp(y0 + 1, 0, m_beardH - 1);

    double tx = x - x0;
    double ty = y - y0;

    auto lerp = [](double a, double b, double t){ return a + (b - a) * t; };

    double c00 = m_beardDensityLUT(y0,x0);
    double c10 = m_beardDensityLUT(y0,x1);
    double c01 = m_beardDensityLUT(y1,x0);
    double c11 = m_beardDensityLUT(y1,x1);

    double v0 = lerp(c00, c10, tx);
    double v1 = lerp(c01, c11, tx);
    double w  = lerp(v0,  v1,  ty);

    return std::clamp(w, 0.0, 1.0);
}

void Renderer::loadTex(const std::string fileName) {

    Image img;  // サイズは loadImage の中で決まる

    if (!img.loadImage(fileName)) {
        std::cerr << "[Renderer::loadTex] failed to load " << fileName << std::endl;
        m_beardW = m_beardH = 0;
        m_beardDensityLUT.resize(0, 0);
        return;
    }

    m_beardW = img.width;
    m_beardH = img.height;
    m_beardDensityLUT.resize(m_beardH, m_beardW);

    for (int y = 0; y < m_beardH; ++y) {
        for (int x = 0; x < m_beardW; ++x) {
            const Eigen::Vector3d& c = img.pixels[y * img.width + x];
            double g = c[0];  // Rチャンネルを密度として使う
            m_beardDensityLUT(y, x) = std::clamp(g, 0.0, 1.0);
        }
    }
}

// ======================= Renderer.cpp 差し替えブロック =======================
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





void  Renderer::Kensyou(const unsigned int &samples) const {
    Image image(camera.getFilm().resolution.x(), camera.getFilm().resolution.y());
    int p_x=image.width/2;
    int p_y=image.height/2;
    Ray ray;
    RayHit hit;
    camera.filmView(p_x, p_y, ray);

    if (hitScene(ray, hit)) {
        const Eigen::Vector3d xi = hit.point;
        runMultilayerTest();

    }



}
void runMultilayerTest()
{
    // 例: 2 層の「表皮 + 真皮」スキン
    std::vector<SSSLayer> skin = {
        // 表皮（薄い・吸収強め）

        SSSLayer{
            Eigen::Vector3d(0.44, 0.22, 0.14),
            Eigen::Vector3d(1.2, 1.3, 1.4),
            0.8,
            1.3,
            0.2
        },


        // 真皮（厚い・散乱強め）

        SSSLayer{
            Eigen::Vector3d(0.0011, 0.0024, 0.014),
            Eigen::Vector3d(0.74, 0.88, 1.01),
            0,
            1.3,
            2
        }

    };

    RWSettings rw;
    rw.numPhotons = 500000; // 時間に応じて増減
    rw.numBins    = 100;
    rw.rMax       = 20.0;

    compareMultilayerMultipoleWithRW(skin, rw, "multilayer_test");


}


void Renderer::computeLocalFrame(const Eigen::Vector3d &w, Eigen::Vector3d &u, Eigen::Vector3d &v) {
    if(fabs(w.x()) > 1e-3)
        u = Eigen::Vector3d::UnitY().cross(w).normalized();
    else
        u = Eigen::Vector3d::UnitX().cross(w).normalized();

    v = w.cross(u);
}
