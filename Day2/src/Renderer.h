//
// Created by kango on 2023/04/03.
//

#ifndef DAY_3_RENDERER_H
#define DAY_3_RENDERER_H


#include <vector>
#include "Body.h"
#include "Camera.h"
#include <random>

class Renderer {
public:
    std::vector<Body> bodies;

    Camera camera;
    Color bgColor;

    /// 乱数生成器
    mutable std::mt19937_64 engine;
    mutable std::uniform_real_distribution<> dist;

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

    void diffuseSample(const Eigen::Vector3d &incidentPoint, const Eigen::Vector3d &normal, Ray &out_Ray) const;

    void Renderer::SSSSample(const Eigen::Vector3d &incidentPoint, const double &radius,const Eigen::Vector3d &normal, Ray &out_Ray) const ;


    static void computeLocalFrame(const Eigen::Vector3d &w, Eigen::Vector3d &u, Eigen::Vector3d &v);
};


#endif //DAY_3_RENDERER_H
