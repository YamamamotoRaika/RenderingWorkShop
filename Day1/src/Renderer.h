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

    Renderer(const std::vector<Body> &bodies, Camera camera, Color bgColor=Color::Zero());

    bool hitScene(const Ray &ray, RayHit &hit) const;

    Image render() const;
};


#endif //DAY_3_RENDERER_H
