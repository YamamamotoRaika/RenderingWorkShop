//
// Created by kango on 2023/04/03.
//

#include "Renderer.h"
#include <iostream>

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
                color=bodies[hit.idx].material.color;
                //std::cout<<hit.idx<<std::endl<<" "<<bodies[hit.idx].show()<<" "<<std::endl<<std::endl;
            }else {
                color=bgColor;
            }


            image.pixels[p_idx] = color;
        }
    }

    return image;
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

//radiusの値って具体的になんだろう？
void Renderer::SSSSample(const Eigen::Vector3d &incidentPoint, const double &radius,const Eigen::Vector3d &normal, Ray &out_Ray) const {
    const double phi = 2.0 * EIGEN_PI * rand();
    const double theta = asin(sqrt(rand()));
    const double _r=radius*rand();


    /// normalの方向をy軸とした正規直交基底を作る
    Eigen::Vector3d u, v;
    computeLocalFrame(normal, u, v);

    const double _x = sin(theta) * cos(phi);
    const double _y = cos(theta);
    const double _z = sin(theta) * sin(phi);

    //incident_pointから距離_r離した点を考える
    //角度は今回Φを流用する。修正筆頭ポイント
    const double d_x =_r*cos(phi);
    const double d_y =_r*sin(phi);

    out_Ray.org = incidentPoint+d_x*u+d_y*v;
}

double A(double eta) {
    double Fdr = -1.440 / (eta * eta) + 0.710 / eta + 0.668 + 0.0636 * eta;
    return (1 + Fdr) / (1 - Fdr);
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
                        Ray sss_ray; RayHit _hit;
                        SSSSample(hit.point, 5.0, hit.normal, sss_ray); // 入射点サンプル

                        // 🔽 スレッドローカル蓄積
                        Color irradiance = Color::Zero();

                        const int inner_samples = 100;
                        std::vector<Color> buffer_irradiance(omp_get_max_threads(), Color::Zero());

#pragma omp parallel
                        {
                            int tid = omp_get_thread_num();
                            for (int j = 0; j < inner_samples / omp_get_num_threads(); ++j) {
                                Ray local_ray;
                                diffuseSample(sss_ray.org, hit.normal, local_ray);

                                RayHit local_hit;
                                if (hitScene(local_ray, local_hit) && bodies[local_hit.idx].isLight()) {
                                    Color contrib = bodies[hit.idx].getKd().cwiseProduct(bodies[local_hit.idx].getEmission());
                                    buffer_irradiance[tid] += contrib;
                                }
                            }
                        }

                        for (const auto &c : buffer_irradiance)
                            irradiance += c;

                        irradiance /= static_cast<double>(inner_samples);

                        reflectRadiance += evaluateBSSRDF(sss_ray.org, hit.point, bodies[hit.idx].material).cwiseProduct(irradiance);
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



void Renderer::computeLocalFrame(const Eigen::Vector3d &w, Eigen::Vector3d &u, Eigen::Vector3d &v) {
    if(fabs(w.x()) > 1e-3)
        u = Eigen::Vector3d::UnitY().cross(w).normalized();
    else
        u = Eigen::Vector3d::UnitX().cross(w).normalized();

    v = w.cross(u);
}
