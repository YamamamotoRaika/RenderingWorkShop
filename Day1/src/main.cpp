#include <iostream>
#include "Body.h"
#include "Camera.h"
#include "Renderer.h"

void intersectTest() {
    const Sphere sphere(1, Eigen::Vector3d::Zero());
    Ray ray(Eigen::Vector3d(0, 0, 10), Eigen::Vector3d(0, 0, -1));
    RayHit hit;
    sphere.hit(ray, hit);

    std::cout << "t:\t" << hit.t << std::endl;
    std::cout << "normal:\t(" << hit.normal.transpose() << ")" << std::endl;
}

void sample() {
    const std::vector<Body> bodies =  {
            Body(Sphere(1.0, Eigen::Vector3d::Zero()), Material(Color(1, 0, 0))),
            Body(Sphere(1.0, Eigen::Vector3d(0, 3, 0)), Material(Color(0, 1, 0))),
            Body(Sphere(1.0, Eigen::Vector3d(0, -3, 0)), Material(Color(0, 0, 1))),
    };

    const Eigen::Vector3d campos(0, 10, 100);
    const Eigen::Vector3d camdir = Eigen::Vector3d(0, 0, 0) - campos;

    const Camera camera(campos, camdir, 320, 9.0 / 16.0, 5);

    /// 背景色はわかりやすく灰色
    const Renderer renderer(bodies, camera, Color(0.1, 0.1, 0.1));
    const auto image = renderer.render().apply_reinhard_extended_tone_mapping().apply_gamma_correction();

    image.save("sample_image.png");
}

int main() {
    std::cout << "Hello, World!" << std::endl;
    intersectTest();
    sample();
    return 0;
}
