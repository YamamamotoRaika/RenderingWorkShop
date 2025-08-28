//
// Created by nakat on 2025/07/09.
//
#include "Triangle.h"

Triangle::Triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c)
    : v0(a), v1(b), v2(c) , n((b-a).cross(c-a).normalized()){}

Triangle::Triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c, const Eigen::Vector3d &n)
    : v0(a), v1(b), v2(c) , n(n){}

bool Triangle::hit(const Ray &ray, RayHit &hit) const {
    const std::array<Eigen::Vector3d, 3> points = {v0, v1, v2};
    const double EPSILON = 1e-8;
    //pawapoを見て
    const Eigen::Vector3d e1=v1-v0;
    const Eigen::Vector3d e2=v2-v0;

    const Eigen::Vector3d alpha = ray.dir.cross(e2);
    double det = e1.dot(alpha);
    //レイが面に対して平衡に入射したとき、falseを返す
    if (-EPSILON  < det && det < EPSILON ) {
        return false;
    }
    double invDet = 1.0f / det;
    Eigen ::Vector3d r=ray.org-v0;
    //ここまでで計算終了、判定を行う

    //0<=u<=1であるかどうか
    const double u=alpha.dot(r)*invDet;
    if (u<0.0 || u>1.0) return false;
    //0<=v<=1かつu+v<=1であるかどうか
    const Eigen::Vector3d beta=r.cross(e1);
    const double v=ray.dir.dot(beta)*invDet;
    if (v<0.0 || v>1.0-u) return false;
    //t>=0であるか(閾値以下であるかも考える)
    const double t=e2.dot(beta)*invDet;
    if (t<1e-6) return false;

    //ここまで行けば交差が確認できるのでRayhit hitに情報を確認したのちtrueを返す。

    hit.t=t;
    hit.point=ray.at(hit.t);
    hit.normal=n;
    if (hit.normal.dot(ray.dir) > 0.0) {
        hit.normal = -hit.normal;
    }



    return true;
}

