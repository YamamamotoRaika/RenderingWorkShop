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
    /// bodiesに光源を追加
    const std::vector<Body> bodies =  {
            Body(Sphere(1.0, Eigen::Vector3d::Zero()), Material(Color(1, 0.1, 0.1), 0.8)),
            Body(Sphere(1.0, Eigen::Vector3d(0, 3, 0)), Material(Color(0.1, 1, 0.1), 0.8)),
            Body(Sphere(1.0, Eigen::Vector3d(0, -3, 0)), Material(Color(0.1, 0.1, 1), 0.8)),
            Body(Sphere(2.0, Eigen::Vector3d(0, 10, 10)), Material(Color(1, 1, 1), 0.8, 10)),
    };

    const Eigen::Vector3d campos(0, 10, 100);
    const Eigen::Vector3d camdir = Eigen::Vector3d(0, 0, 0) - campos;

    const Camera camera(campos, camdir, 320, 9.0 / 16.0, 5);

    /// 背景色はわかりやすく灰色
    const Renderer renderer(bodies, camera, Color(0.1, 0.1, 0.1));
    const unsigned int samples = 100;
    const auto image1 = renderer.render();
    const auto image2 = renderer.directIlluminationRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();

    image1.save("sample_image.png");
    image2.save("sample_image2.png");
}

void roomRenderingSample() {
    const auto room_r = 1e5;
    const auto floor_color = Color(0.75, 0.75, 0.75);
    const std::vector<Body> room_walls {
            Body(Sphere(room_r, (room_r - 50) * Eigen::Vector3d::UnitX()), Material(Color(0.25, 0.25, 0.75), 0.8)),
            Body(Sphere(room_r, -(room_r - 50) * Eigen::Vector3d::UnitX()), Material(Color(0.25, 0.75, 0.25), 0.8)),
            Body(Sphere(room_r, (room_r - 30) * Eigen::Vector3d::UnitY()), Material(floor_color, 0.8)),
            Body(Sphere(room_r, -(room_r - 35) * Eigen::Vector3d::UnitY()), Material(floor_color, 0.8)),
            Body(Sphere(room_r, (room_r - 30) * Eigen::Vector3d::UnitZ()), Material(floor_color, 0.8)),
    };

    std::vector<Body> bodies {
            Body(Sphere(16.5, Eigen::Vector3d(25, -14.5, -10)), Material(Color(0.75, 0.25, 0.25), 0.8)),
            Body(Sphere(16.5, Eigen::Vector3d(-23, -14.5, 15)), Material(Color(0.99, 0.99, 0.99), 0.8)),
    };

    const std::vector<Body> lights {
            Body(Sphere(5, Eigen::Vector3d(0, 28, 0)), Material(Color(1, 1, 1), 1.0, 10))
    };

    for(const auto & room_wall : room_walls) {
        bodies.push_back(room_wall);
    }

    for(const auto & light : lights) {
        bodies.push_back(light);
    }

    const Eigen::Vector3d campos(0, 0, 80);
    const Eigen::Vector3d camdir = Eigen::Vector3d(0, 0, 0) - campos;

    const Camera camera(campos, camdir, 540, 4.0 / 3.0, 60, 45);

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
