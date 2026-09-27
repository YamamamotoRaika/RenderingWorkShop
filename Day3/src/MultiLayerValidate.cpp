//
// Created by nakat on 2025/11/19.
//
// MultiLayerValidate.cpp
#include <Eigen/Core>
#include <vector>
#include <random>
#include <fstream>
#include <iostream>
#include <cmath>

// あなたの既存のヘッダ
#include "MultiLayerBSSRDF.h"
#include "MultiLayerValidate.h"

#include <iomanip>
// SSSLayer, Rd_multilayer_multipole を定義しているヘッダ
// もしここに無ければ、SSSLayer の宣言だけこのファイルにもってきてください。

// --------------------------------------
// RW 用の設定
// --------------------------------------


// 半径 r に対応するビン index を返す（範囲外なら -1）
inline int radialBin(double r, double rMax, unsigned int numBins) {
    if (r < 0.0 || r >= rMax) return -1;
    const double dr = rMax / numBins;
    int idx = static_cast<int>(r / dr);
    if (idx >= static_cast<int>(numBins)) idx = static_cast<int>(numBins) - 1;
    return idx;
}

// 方向を一様球面分布でサンプル
inline Eigen::Vector3d sampleUniformSphere(std::mt19937_64 &rng) {
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    double z   = 2.0 * uni(rng) - 1.0;
    double phi = 2.0 * EIGEN_PI * uni(rng);
    double r   = std::sqrt(std::max(0.0, 1.0 - z * z));
    return Eigen::Vector3d(r * std::cos(phi), r * std::sin(phi), z);
}// cosThetaI : 入射方向と法線のなす角の cos 値（表側から見て 0〜1）
// etaI      : 入射側の屈折率 (n_in)
// etaT      : 透過側の屈折率 (n_out)
//
// 戻り値    : 境界での反射率 F（確率として 0〜1）
inline double fresnelDielectric(double cosThetaI, double etaI, double etaT)
{
    // 負の値が来ても大丈夫なように
    cosThetaI = std::clamp(cosThetaI, -1.0, 1.0);

    // 内外の入れ替え（法線の向きと cos の符号が反転している場合）
    bool entering = cosThetaI > 0.0;
    if (!entering) {
        std::swap(etaI, etaT);
        cosThetaI = std::fabs(cosThetaI);
    }

    // Snell の法則から sin^2(theta_t) を計算
    double eta  = etaI / etaT;
    double sin2ThetaI = std::max(0.0, 1.0 - cosThetaI * cosThetaI);
    double sin2ThetaT = eta * eta * sin2ThetaI;

    // 全反射
    if (sin2ThetaT >= 1.0) {
        return 1.0;
    }

    double cosThetaT = std::sqrt(std::max(0.0, 1.0 - sin2ThetaT));

    // フレネル反射係数 (Rs, Rp) を計算して平均
    double Rs_num = etaI * cosThetaI - etaT * cosThetaT;
    double Rs_den = etaI * cosThetaI + etaT * cosThetaT;
    double Rs = (Rs_num / Rs_den) * (Rs_num / Rs_den);

    double Rp_num = etaT * cosThetaI - etaI * cosThetaT;
    double Rp_den = etaT * cosThetaI + etaI * cosThetaT;
    double Rp = (Rp_num / Rp_den) * (Rp_num / Rp_den);

    return 0.5 * (Rs + Rp);
}


// 単純化した Henyey-Greenstein → g=0 のとき一様球面
inline Eigen::Vector3d sampleHG(const Eigen::Vector3d &wo, double g, std::mt19937_64 &rng) {
    if (std::abs(g) < 1e-6) {
        return sampleUniformSphere(rng);
    }
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    double u1 = uni(rng);
    double u2 = uni(rng);

    double cosTheta;
    if (std::abs(g) < 1e-3) {
        cosTheta = 1.0 - 2.0 * u1;
    } else {
        double sq = (1.0 - g * g) / (1.0 - g + 2.0 * g * u1);
        cosTheta = (1.0 + g * g - sq * sq) / (2.0 * g);
    }
    double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    double phi = 2.0 * EIGEN_PI * u2;

    // ローカル座標
    Eigen::Vector3d w(0, 0, 1);
    Eigen::Vector3d u, v;
    if (std::fabs(wo.z()) > 0.1) {
        u = wo.cross(Eigen::Vector3d::UnitX()).normalized();
    } else {
        u = wo.cross(Eigen::Vector3d::UnitZ()).normalized();
    }
    v = wo.cross(u);

    // ローカルからワールドへ
    Eigen::Vector3d wi = (std::cos(phi) * sinTheta) * u +
                         (std::sin(phi) * sinTheta) * v +
                         (cosTheta) * wo;
    return wi.normalized();
}

// z から属する layer index を返す（z>=0: surface, z>0 が内部）
inline int findLayer(double z, const std::vector<SSSLayer> &layers, double totalThickness) {
    if (z < 0.0 || z >= totalThickness) return -1;
    double current = 0.0;
    for (int i = 0; i < static_cast<int>(layers.size()); ++i) {
        double next = current + layers[i].thickness;
        if (z >= current && z < next) return i;
        current = next;
    }
    return -1;
}

// --------------------------------------
// ランダムウォーク 1 本を追跡
//   戻り値: 成功して表面から出たかどうか
// --------------------------------------
bool randomWalkMultiLayer(const std::vector<SSSLayer> &layers,
                          const RWSettings &rwOpt,
                          double totalThickness,
                          std::mt19937_64 &rng,
                          double &radius,
                          Eigen::Vector3d &throughput)
{
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    // 初期状態：原点直下から内部に少し入った位置、下向き
    Eigen::Vector3d p(0.0, 0.0, rwOpt.zStartOffset);
    Eigen::Vector3d w(0.0, 0.0, 1.0); // 下向き

    throughput = Eigen::Vector3d::Ones();

    const int maxBounces = 10000;

    auto layerZRange = [&](int li, double &zMin, double &zMax) {
        zMin = 0.0;
        for (int i = 0; i < li; ++i) {
            zMin += layers[i].thickness;
        }
        zMax = zMin + layers[li].thickness;

    };

    for (int bounce = 0; bounce < maxBounces; ++bounce) {

        int li = findLayer(p.z(), layers, totalThickness);
        if (li < 0) {
            // 外に出ている
            return false;
        }
        const auto &L = layers[li];

        Eigen::Vector3d sigma_s = L.sigma_s;
        Eigen::Vector3d sigma_a = L.sigma_a;
        Eigen::Vector3d sigma_t = sigma_s + sigma_a;

        // 簡単のため 3 チャンネルの平均で free-flight をサンプル
        double sigma_t_bar = (sigma_t[0] + sigma_t[1] + sigma_t[2]) * (1.0 / 3.0);
        if (sigma_t_bar <= 0.0) return false;

        double u = uni(rng);
        double s = -std::log(std::max(1e-12, 1.0 - u)) / sigma_t_bar;

        // 現在層の z 範囲 [zMin, zMax]
        double zMin, zMax;
        layerZRange(li, zMin, zMax);

        // この層の上下境界までの距離
        double tToTopLayer    = std::numeric_limits<double>::infinity();
        double tToBottomLayer = std::numeric_limits<double>::infinity();

        if (w.z() < -1e-12) {
            tToTopLayer = (zMin - p.z()) / w.z();   // zMin <= p.z(), w.z()<0 → t>0
        }
        if (w.z() >  1e-12) {
            tToBottomLayer = (zMax - p.z()) / w.z(); // zMax >= p.z(), w.z()>0 → t>0
        }

        double tHit = std::min(tToTopLayer, tToBottomLayer);
        bool hitTopBoundary    = (tToTopLayer    < tToBottomLayer);
        bool hitBottomBoundary = !hitTopBoundary;

        if (s >= tHit) {
            // --- 境界に到達してしまう場合 → 境界まで進める ---
            double travel = tHit;
            p += travel * w;

            // Beer-Lambert（境界までの吸収）
            Eigen::Vector3d Tr;
            Tr[0] = std::exp(-sigma_a[0] * travel);
            Tr[1] = std::exp(-sigma_a[1] * travel);
            Tr[2] = std::exp(-sigma_a[2] * travel);
            throughput = throughput.cwiseProduct(Tr);

            // どの境界か（外界 or 層間）を判定
            bool isTopSurface    = hitTopBoundary    && (li == 0);
            bool isBottomSurface = hitBottomBoundary && (li == (int)layers.size() - 1);

            // 境界の法線と屈折率
            Eigen::Vector3d n;   // 境界「外側」（行き先側）を向く法線
            double eta_in  = L.eta;
            double eta_out = 1.0;

            if (hitTopBoundary) {
                // 上側へ
                n = Eigen::Vector3d(0, 0, -1);
                if (!isTopSurface) {
                    // 層間境界 → 上の層の屈折率
                    eta_out = layers[li - 1].eta;
                }
            } else {
                // 下側へ
                n = Eigen::Vector3d(0, 0, 1);
                if (!isBottomSurface) {
                    // 層間境界 → 下の層の屈折率
                    eta_out = layers[li + 1].eta;
                }
            }

            // 入射方向は「境界に向かう」向きなので -w を使う
            double cos_i = (-w).dot(n);
            // Fresnel（関数の仕様によるが、PBRT 風を想定）
            double F = fresnelDielectric(cos_i, eta_in, eta_out);
            F = std::clamp(F, 0.0, 1.0);

            double uF = uni(rng);
            if (uF < F) {
                // --- 反射：層は変わらず、方向だけ鏡面反射 ---
                w = (w - 2.0 * w.dot(n) * n).normalized();
                p += 1e-4 * w; // 境界から少し離す
                continue;      // 次の bounce へ（新しい free-flight をサンプル）
            } else {
                // --- 透過 ---
                if (isTopSurface) {
                    // 上面から外へ → サンプル確定
                    // （w は今のままでも OK：出射方向に依存する項は
                    //    後で cosθ をかける or 既に考慮済みなら不要）
                    radius = std::sqrt(p.x() * p.x() + p.y() * p.y());
                    return true;
                }
                if (isBottomSurface) {
                    // 下面から外へ
                    //   ・黒バック板 → 吸収終了
                    //   ・下に別の媒質（例えば骨など）を入れたいなら
                    //     ここでその媒質に入って続行、という拡張も可能
                    return false;
                }

                // --- 層間境界を透過 → 隣の層へ ---
                // ここで本当は Snell の法則で屈折方向を計算するのが一番厳密ですが、
                // 多重散乱が支配的な BSSRDF では屈折による曲がりは 2 次効果なので、
                // まずは方向をそのままにして「層だけ切り替える」近似でも十分です。
                //
                // ※もし屈折も入れたいなら、ここで refract(...) を呼んで w を更新する。

                p += 1e-4 * w; // 次の層側に少しだけ押し出す
                // 次の bounce で findLayer(p.z(), ...) し直されるので、
                // 層インデックスはそこで自動的に更新される。
                continue;
            }
        } else {
            // --- 媒質内で散乱 ---
            p += s * w;

            // Beer-Lambert（吸収）
            Eigen::Vector3d Tr;
            Tr[0] = std::exp(-sigma_a[0] * s);
            Tr[1] = std::exp(-sigma_a[1] * s);
            Tr[2] = std::exp(-sigma_a[2] * s);

            Eigen::Vector3d albedo = sigma_s.cwiseQuotient(sigma_t);
            throughput = throughput.cwiseProduct(Tr.cwiseProduct(albedo));

            // ロシアンルーレット
            double maxComp = throughput.maxCoeff();
            if (maxComp < 1e-3) {
                double q = std::min(0.9, maxComp * 10.0);
                if (uni(rng) > q) return false;
                throughput /= q;
            }

            // 散乱方向をサンプル
            w = sampleHG(w, L.g, rng);
            w = sampleHG(w,0, rng);
        }
    }

    // バウンス上限
    return false;
}

// --------------------------------------
// マルチポール vs RW の比較本体
// --------------------------------------
void compareMultilayerMultipoleWithRW(const std::vector<SSSLayer> &layers,
                                      const RWSettings &rwOpt,
                                      const std::string &csvPrefix)
{
    const double rMax = rwOpt.rMax;
    const unsigned int numBins = rwOpt.numBins;
    const double dr = rMax / numBins;

    // 全厚さ
    double totalThickness = 0.0;
    for (const auto &L : layers) {

        totalThickness += L.thickness*2;
    }
    std::cout<<totalThickness<<std::endl;
    // ヒストグラム
    std::vector<Eigen::Vector3d> fluxRW(numBins, Eigen::Vector3d::Zero());

    // ---- ランダムウォーク ----
    std::mt19937_64 rng(42);

    for (unsigned int i = 0; i < rwOpt.numPhotons; ++i) {
        double r;
        Eigen::Vector3d T;

        std::vector<SSSLayer>Layers=layers;
        for (SSSLayer &L : Layers) {
            L.setThickness(L.thickness*2);

        }

        if (randomWalkMultiLayer(Layers, rwOpt, totalThickness, rng, r, T)) {
            int idx = radialBin(r, rMax, numBins);
            if (idx >= 0) {
                fluxRW[idx] += T;
            }
        }
    }

    // Rd_RW(r) へ変換
    std::vector<Eigen::Vector3d> RdRW(numBins, Eigen::Vector3d::Zero());
    for (unsigned int i = 0; i < numBins; ++i) {
        double r0  = i * dr;
        double r1  = (i + 1) * dr;
        double rmid = 0.5 * (r0 + r1);
        double areaRing = 2.0 * EIGEN_PI * rmid * dr;

        if (areaRing > 0.0) {
            RdRW[i] = fluxRW[i] / (static_cast<double>(rwOpt.numPhotons) * areaRing);
        }
    }

    // ---- マルチポールで Rd(r) を計算 ----
    std::vector<Eigen::Vector3d> RdMP(numBins, Eigen::Vector3d::Zero());

    std::vector<SSSLayer> xxxx = layers;
    for (auto& xx : xxxx) xx.setDerived();


    for (unsigned int i = 0; i < numBins; ++i) {
        double r0   = i * dr;
        double r1   = (i + 1) * dr;
        double rmid = 0.5 * (r0 + r1);
        MLProfileSettings opt;
        opt.imagesPerSide = 8;

        double C_boundry=(1-xxxx[0].Fdr);



        RdMP[i] = C_boundry*C_boundry*Rd_multipole_multilayer(layers,rmid,opt);


    }

    // ---- p_A(r) 用の正規化係数 Z を計算 ----
    auto integrateEnergy = [&](const std::vector<Eigen::Vector3d> &Rd) {
        Eigen::Vector3d Z = Eigen::Vector3d::Zero();
        for (unsigned int i = 0; i < numBins; ++i) {
            double r0   = i * dr;
            double r1   = (i + 1) * dr;
            double rmid = 0.5 * (r0 + r1);
            double areaRing = 2.0 * EIGEN_PI * rmid * dr;
            Z += Rd[i] * areaRing;
        }
        return Z;
    };

    Eigen::Vector3d Zrw = integrateEnergy(RdRW);
    Eigen::Vector3d Zmp = integrateEnergy(RdMP);

    // 1) Rd プロファイルを書き出し
    {
        std::ofstream ofs(csvPrefix + "_Rd.csv");
        ofs << "# r, Rd_rw_R, Rd_rw_G, Rd_rw_B, Rd_mp_R, Rd_mp_G, Rd_mp_B\n";
        for (unsigned int i = 0; i < numBins; ++i) {
            double r0   = i * dr;
            double r1   = (i + 1) * dr;
            double rmid = 0.5 * (r0 + r1);
            ofs << rmid << ","
                << RdRW[i][0] << "," << RdRW[i][1] << "," << RdRW[i][2] << ","
                << RdMP[i][0] << "," << RdMP[i][1] << "," << RdMP[i][2] << "\n";
        }
        std::cout << "Saved Rd profile to " << csvPrefix << "_Rd.csv\n";
    }

    // 2) p_A(r)（エネルギー正規化された半径方向 pdf）を書き出し
    {
        std::ofstream ofs(csvPrefix + "_pA.csv");
        ofs << "# r, pA_rw_R, pA_rw_G, pA_rw_B, pA_mp_R, pA_mp_G, pA_mp_B\n";
        for (unsigned int i = 0; i < numBins; ++i) {
            double r0   = i * dr;
            double r1   = (i + 1) * dr;
            double rmid = 0.5 * (r0 + r1);
            double areaRing = 2.0 * EIGEN_PI * rmid * dr;

            Eigen::Vector3d p_rw = Eigen::Vector3d::Zero();
            Eigen::Vector3d p_mp = Eigen::Vector3d::Zero();
            if (areaRing > 0.0) {
                if (Zrw[0] > 0.0) p_rw[0] = RdRW[i][0] * areaRing / Zrw[0];
                if (Zrw[1] > 0.0) p_rw[1] = RdRW[i][1] * areaRing / Zrw[1];
                if (Zrw[2] > 0.0) p_rw[2] = RdRW[i][2] * areaRing / Zrw[2];

                if (Zmp[0] > 0.0) p_mp[0] = RdMP[i][0] * areaRing / Zmp[0];
                if (Zmp[1] > 0.0) p_mp[1] = RdMP[i][1] * areaRing / Zmp[1];
                if (Zmp[2] > 0.0) p_mp[2] = RdMP[i][2] * areaRing / Zmp[2];
            }

            ofs << rmid << ","
                << p_rw[0] << "," << p_rw[1] << "," << p_rw[2] << ","
                << p_mp[0] << "," << p_mp[1] << "," << p_mp[2] << "\n";
        }
        std::cout << "Saved p_A(r) to " << csvPrefix << "_pA.csv\n";
    }


}
