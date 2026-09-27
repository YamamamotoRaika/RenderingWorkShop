//
// Created by nakat on 2025/11/19.
//

#ifndef MULTILAYERVALIDATE_H
#define MULTILAYERVALIDATE_H
// MultiLayerValidate.h

struct RWSettings {
    unsigned int numPhotons   = 200000; // サンプル数（時間と相談して調整）
    unsigned int numBins      = 100;    // r の分割数
    double       rMax         = 20.0;   // 最大半径 [mm]
    double       zStartOffset = 1e-4;   // 表面から少しだけ内部から開始
};
#pragma once

void runMultilayerTest();
void compareMultilayerMultipoleWithRW(const std::vector<SSSLayer> &layers,
                                      const RWSettings &rwOpt,
                                      const std::string &csvPrefix);

#endif //MULTILAYERVALIDATE_H
