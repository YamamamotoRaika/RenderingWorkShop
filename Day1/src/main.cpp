#include <iostream>
#include "Body.h"
#include "Camera.h"
#include "Renderer.h"

/// \衝突判定のテスト
/// 球を原点に用意してレイを飛ばします
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
            Body(Sphere(3.0, Eigen::Vector3d(0, 0, -10)), Material(Color(1, 1, 1))),
    };

    const Eigen::Vector3d campos(0, 0, 100);
    const Eigen::Vector3d camdir = Eigen::Vector3d(0, 0, 0) - campos;

    const Camera camera(campos, camdir, 320, 9.0 / 16.0, 5);

    /// 背景色はわかりやすく灰色
    const Renderer renderer(bodies, camera, Color(0.5, 0.5, 0.5));
    const auto image = renderer.render();

    image.save("day1_goal.png");
}

void roomRenderingSample() {
    const auto room_r = 1e5;
    const auto floor_color = codeToColor("#f9c89b");
    const std::vector<Body> room_walls {
            Body(Sphere(room_r, (room_r - 50) * Eigen::Vector3d::UnitX()), Material(codeToColor("#2f5d50"))),
            Body(Sphere(room_r, -(room_r - 50) * Eigen::Vector3d::UnitX()), Material(codeToColor("#00a3af"))),
            Body(Sphere(room_r, (room_r - 30) * Eigen::Vector3d::UnitY()), Material(floor_color)),
            Body(Sphere(room_r, -(room_r - 40) * Eigen::Vector3d::UnitY()), Material(floor_color)),
            Body(Sphere(room_r, (room_r - 30) * Eigen::Vector3d::UnitZ()), Material(floor_color)),
    };

    std::vector<Body> bodies {
            Body(Sphere(16.5, Eigen::Vector3d(25, -14.5, -10)), Material(Color(0.75, 0.25, 0.25))),
            Body(Sphere(16.5, Eigen::Vector3d(-23, -14.5, 7)), Material(Color(0.99, 0.99, 0.99))),
    };

    const std::vector<Body> lights {
            Body(Sphere(5, Eigen::Vector3d(0, 34.8, 0)), Material(codeToColor("#e597b2")))
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
    roomRenderingSample();
    sample();
    return 0;
}
