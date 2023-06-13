//
// Created by kango on 2023/04/14.
//

#ifndef DAY_3_MATERIAL_H
#define DAY_3_MATERIAL_H

#include <utility>

#include "Image.h"

struct Material {
    Color color;

public:
    Material(Color color) : color(std::move(color)) {}
};

#endif //DAY_3_MATERIAL_H
