//
// Created by nakat on 2025/07/09.
//
#include "Triangle.h"

Triangle::Triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c)
    : v0(a), v1(b), v2(c) , n((b-a).cross(c-a).normalized()){}

Triangle::Triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c, const Material &m)
    : v0(a), v1(b), v2(c) , material(m),n((b-a).cross(c-a).normalized()){}

Triangle::Triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c, const Eigen::Vector3d &n)
    : v0(a), v1(b), v2(c) , n(n){}

Triangle::Triangle(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c,
    const Eigen::Vector2d &d, const Eigen::Vector2d &e, const Eigen::Vector2d &f,const Eigen::Vector3d &n)
    : v0(a), v1(b), v2(c) ,uv0(d), uv1(e), uv2(f) ,n(n){}

Eigen::Vector3d Triangle::centroid() const {
    return (v0 + v1 + v2) / 3.0;
}
AABB Triangle::computeAABB() const {
    Eigen::Vector3d minPt = v0.cwiseMin(v1).cwiseMin(v2);
    Eigen::Vector3d maxPt = v0.cwiseMax(v1).cwiseMax(v2);
    return AABB(minPt, maxPt);
}
Eigen::Vector2d Triangle::getUV( const Eigen::Vector3d point) const {

    // v0 を原点とした座標系に変換
    Eigen::Vector3d v0v1 = v1 - v0;
    Eigen::Vector3d v0v2 =v2 - v0;
    Eigen::Vector3d v0p  = point  - v0;

    // バリセントリック (b1,b2) を解く
    double d00 = v0v1.dot(v0v1);
    double d01 = v0v1.dot(v0v2);
    double d11 = v0v2.dot(v0v2);
    double d20 = v0p .dot(v0v1);
    double d21 = v0p .dot(v0v2);

    double denom = d00 * d11 - d01 * d01;
    if (std::abs(denom) < 1e-12) {
        // ほぼ退化三角形のとき、とりあえず v0 の uv を返す
        return uv0;
    }

    double b1 = (d11 * d20 - d01 * d21) / denom; // 頂点1の重み
    double b2 = (d00 * d21 - d01 * d20) / denom; // 頂点2の重み
    double b0 = 1.0 - b1 - b2;                  // 頂点0の重み

    // 同じ重みで UV を補間
    Eigen::Vector2d uv =
        b0 * uv0 +
        b1 * uv1 +
        b2 * uv2;

    return uv;


}

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
    if (t<1e-4) return false;

    //ここまで行けば交差が確認できるのでRayhit hitに情報を確認したのちtrueを返す。

    hit.t=t;
    hit.point=ray.at(hit.t);
    hit.normal=n;
    if (hit.normal.dot(ray.dir) > 0.0) {
        hit.normal = -hit.normal;
    }



    return true;
}
double Triangle::area() const {
    return 0.5 * (v1 - v0).cross(v2 - v0).norm();
}
Material Triangle::getMaterial() const {
    return  material;
}
Eigen::Vector3d Triangle::getKd() const {
    return material.color*material.kd;
}
Eigen::Vector3d Triangle::getEmission() const {
    return material.color*material.emission;
}
bool Triangle::isLight() const {
    return material.emission > 0.0;
}