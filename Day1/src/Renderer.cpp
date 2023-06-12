//
// Created by kango on 2023/04/03.
//

#include "Renderer.h"

#include <utility>
#include <iostream>

Renderer::Renderer(const std::vector<Body> &bodies, Camera camera, Color bgColor)
: bodies(bodies), camera(std::move(camera)), bgColor(std::move(bgColor)), engine(0), dist(0, 1) {}

double Renderer::rand() const {
    return dist(engine);
}

bool Renderer::hitScene(const Ray &ray, RayHit &hit) const {
    hit.t = DBL_MAX;
    hit.idx = -1;
    for(int i = 0; i < bodies.size(); i++) {
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

void Renderer::directionalSample(const RayHit &hit, Color &out_color) const {
    Ray _ray; RayHit _hit;
    diffuseSample(hit.point, hit.normal, _ray);
    if(hitScene(_ray, _hit)) {
        out_color = bodies[hit.idx].getKd().cwiseProduct(bodies[_hit.idx].getEmission());
    } else {
        out_color = Color::Zero();
    }
}

void Renderer::diffuseSample(const Eigen::Vector3d &incidentPoint, const Eigen::Vector3d &normal, Ray &out_Ray) const {
    const double phi = 2.0 * EIGEN_PI * dist(engine);
    const double theta = asin(sqrt(dist(engine)));

    /// normalの方向をy軸とした正規直交基底を作る
    Eigen::Vector3d u, v;
    computeLocalFrame(normal, u, v);

    const double _x = sin(theta) * cos(phi);
    const double _y = cos(theta);
    const double _z = sin(theta) * sin(phi);

    out_Ray.dir = _x * u + _y * normal + _z * v;
    out_Ray.org = incidentPoint;
}

void Renderer::computeLocalFrame(const Eigen::Vector3d &w, Eigen::Vector3d &u, Eigen::Vector3d &v) {
    if(fabs(w.x()) > 1e-3)
        u = Eigen::Vector3d::UnitY().cross(w).normalized();
    else
        u = Eigen::Vector3d::UnitX().cross(w).normalized();

    v = w.cross(u);
}
