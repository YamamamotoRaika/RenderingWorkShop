//
// Created by kango on 2023/04/03.
//

#ifndef DAY_3_RENDERER_H
#define DAY_3_RENDERER_H


#include <vector>
#include "Body.h"
#include "Camera.h"
#include <random>
#include <unordered_map>
#include <memory>
#include <mutex>

struct BSSRDFRadialCDF { /* ←前回のテーブル型そのまま */
    std::vector<double> rGrid, cdf;
    double rMax = 0.0, Z = 1.0;
    void build(const Material& material, int N = 1024);
    double sampleR(double u) const;
    double radialPdf(double r) const; // p_r(r) = 2π r Rd(r) / Z
    void buildSS(const Material& material, int N = 1024);

};


class Renderer {
public:
    std::vector<Body> bodies;

    Camera camera;
    Color bgColor;

    /// 乱数生成器
    mutable std::mt19937_64 engine;
    mutable std::uniform_real_distribution<> dist;
    mutable std::unordered_map<const Material*, std::shared_ptr<BSSRDFRadialCDF>> sssRadialCache;
    mutable std::unordered_map<const Material*, std::shared_ptr<BSSRDFRadialCDF>> sssRadialCacheMS;
    mutable std::unordered_map<const Material*, std::shared_ptr<BSSRDFRadialCDF>> sssRadialCacheSS;
    mutable std::mutex sssCacheMtx;

    const BSSRDFRadialCDF& getRadialTable(const Material& m) const;
    const BSSRDFRadialCDF& getRadialTableMS(const Material& m) const;
    const BSSRDFRadialCDF& getRadialTableSS(const Material& m) const;


    Renderer(const std::vector<Body> &bodies, Camera camera, Color bgColor=Color::Zero());

    double rand() const;

    bool hitScene(const Ray &ray, RayHit &hit) const;

    Image render() const;




    Image directIlluminationRender(const unsigned int &samples) const;

    Image _directIlluminationRender(const unsigned int &samples) const;
    Image SSSdirectIlluminationRender(const unsigned int &samples) const;
    Image _SSSdirectIlluminationRender(const unsigned int &samples) const;
    Color Renderer::evaluateBSSRDF(const Eigen::Vector3d& xi,
                                const Eigen::Vector3d& xo,
                                const Material& material) const;
    Color Renderer::_evaluateBSSRDF(const Eigen::Vector3d& xi,
                                const Eigen::Vector3d& xo,
                                const Material& material) const;

    static double evaluateBSSRDFScalarAtDistance(double d, const Material& material);
    static double evaluateBSSRDFScalarAtDistanceSS(double d, const Material& material); // SS

    void diffuseSample(const Eigen::Vector3d &incidentPoint, const Eigen::Vector3d &normal, Ray &out_Ray) const;

    void Renderer::SSSSample(const Eigen::Vector3d &incidentPoint, const double &radius, double &r,const Eigen::Vector3d &normal, Ray &out_Ray) const ;
    void Renderer::_SSSSample(const Eigen::Vector3d &incidentPoint,
                         const Eigen::Vector3d &normal,
                         const Material& material,
                         Ray &out_Ray,double &out_pdfA) const ;
    void _SSSSampleMixed(const Eigen::Vector3d& incidentPoint,
                       const Eigen::Vector3d& normal,
                       const Material& material,
                       Ray& out_Ray,
                       double& out_pdfA) const;
    bool Renderer::projectToSurface(const Eigen::Vector3d& xop,
                                const Eigen::Vector3d& n_hint,
                                int bodyId,
                                Eigen::Vector3d& xo,
                                Eigen::Vector3d& n_o) const;
   Image Renderer::KAI_SSSdirectIlluminationRender(const unsigned int &samples) const ;
//検証用コード
    inline double Renderer::urand() const;
    Eigen::Vector3d Renderer::sampleIsotropicDir(double u1, double u2) const;
    double Renderer::sampleFreeFlight(double u, double sigma_t) const;
    bool Renderer::refract(const Eigen::Vector3d &wi, const Eigen::Vector3d &n, double n1, double n2,
                       Eigen::Vector3d &wt) const ;
    double Renderer::fresnelDielectric_R(double n1, double n2, double cosThetaI) const ;
    inline Eigen::Vector3d Renderer::offsetP(const Eigen::Vector3d& p, const Eigen::Vector3d& n) const ;
    Image Renderer::ReferenceSSSRandomWalkRender(const unsigned int &spp) const;

    // 1点入射のラジアル・プロファイルをCSVへ出力（ランダムウォークの経験Rdと理論Rdの比較）
    // - xi, ni : 入射点とその法線（半無限平板の同一点を推奨）
    // - bodyId : xi が載っている Body の id（同一オブジェクトに限定するため必須）
    // - samples: 例) 1'000'000 など
    // - rMax   : 半径上限（例: 0.02 = 2cm。材質/スケールに合わせて）
    // - bins   : ヒストグラムビン数（例: 300）
    // - csvPath: 出力先
    void ProfileRadialBSSRDF_RW(const Eigen::Vector3d& xi,
                                const Eigen::Vector3d& ni,
                                const Material& mat,
                                int bodyId,
                                size_t samples,
                                double rMax,
                                int bins,
                                const std::string& csvPath) const;
    // 1点入射：混合(SS+MS)の理論プロファイルをCSV出力
    void ProfileCompare_Rw_vs_MixedTheory(const Eigen::Vector3d& xi,
                                   const Eigen::Vector3d& ni,
                                   const Material& mat,
                                   int bodyId,
                                   size_t rw_samples,
                                   size_t ss_samples,   // SS用のサンプル数（例: 1e6）
                                   double rMax,
                                   int bins,
                                   const std::string& csvPath) const;
    // 厚さ[m]を与える（半無限にしたければ thickness<=0）
    // M は像のサマンド数（0=ディポール、1..3 で多くの薄板は十分）
    Color evaluateBSSRDF_Multipole(const Eigen::Vector3d& xi,
                                   const Eigen::Vector3d& xo,
                                   const Material& material,
                                   double thickness /*T*/,
                                   double r_max,
                                   int    M /*images on each side*/) const;
    static double evaluateBSSRDF_MultipoleAtDistance(const double r,
                                   const Material& material,
                                   double thickness /*T*/,
                                   double r_max,
                                   int    M /*images on each side*/) ;

      void _ProfileRadialBSSRDF_RW(const Eigen::Vector3d& xi,
                                 const Eigen::Vector3d& ni,
                                 const Material& mat,
                                 int bodyId,
                                 size_t samples,
                                 double thickness,   // [m]
                                 double rMax,        // [m]
                                 int bins,
                                 const std::string& csvPath) const;

    double computeRhoBSSRDF(const Eigen::Vector3d& xi, const Material& mat,
                        double rMax, int N, bool includeFresnelOutside) const;

    void Renderer::Kensyou(const unsigned int &samples) const;




    static void computeLocalFrame(const Eigen::Vector3d &w, Eigen::Vector3d &u, Eigen::Vector3d &v);
};


#endif //DAY_3_RENDERER_H
