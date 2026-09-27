//
// Created by nakat on 2025/10/30.
//

#ifndef MESH_H
#define MESH_H



#include "Eigen/Dense"
#include "Ray.h"
#include  "Triangle.h"

#include  "AABB.h"
#include "Material.h"
#include "Image.h"

// BVHノード構造
struct BVHNode {
    AABB bbox;               // このノードのバウンディングボックス
    BVHNode* left = nullptr; // 左子ノード
    BVHNode* right = nullptr;// 右子ノード
    int start = 0;           // 三角形の開始インデックス（leaf用）
    int range = 0;           // 三角形数（leaf用）
    bool isLeaf = false;

    ~BVHNode() {
        delete left;
        delete right;
    }
};

class Mesh {

private:
    mutable int triId=-1;
    std::vector<AABB> triangleBBoxes;


    BVHNode* buildRecursive(int start, int end, int depth);
    bool hitBVH(const Ray& ray, RayHit& hit, const BVHNode* node) const;
    bool hitBVH(const Ray& ray, RayHit& hit,
                const BVHNode* node, int ignoreTriId) const;

    // デバッグ: ルートAABBのみで判定して三角形は総当り
    AABB meshBBox;
    bool hitSingleBBox(const Ray& ray, RayHit& hit) const;

    void buildBoundsOnly(); // ルートAABBだけ構築

public:
    std::vector<Triangle> triangles;
    std::vector<int>      primIndices;
    std::vector<AABB>     triBBox;     // 各三角形のAABB（同じ順序）
    BVHNode* bvhRoot = nullptr;
    Image m_texture;
    bool m_hasTexture = false;









    Mesh() = default;
    Mesh(const std::vector<Triangle> &triangles);



    bool hit(const Ray& ray, RayHit& hit, int ignoreTriId = -1) const;

    Eigen::Vector3d getKd() const;

    bool isLight() const;
    Material getMaterial() const;
    void getNormal() const;
    Eigen::Vector3d getEmission() const;
    Eigen::Vector2d getUV(const Eigen::Vector3d point) const;

    Eigen::Vector3d getPos() const;
    Eigen::Vector3d getSize() const;

    // スケーリング
    void scale(double s);





    void buildBVH();
    void setTriId(int i) const;
    void setTexture(const std::string& fileName);
    bool loadTexture(const std::string& fileName);
    Color getTexture(const Eigen::Vector2d& uv) const;
    void setBSSRDFParams(double scale,double eta,Eigen::Vector3d sigma_a,Eigen::Vector3d sigma_s);
    void addBeardFromRootsCSV(const std::string& rootsPath,
                              double hairLen,
                              double hairRad,
                              double ySurface,
                              bool rootsAreMeters,
                              const Material& hairMat);
    void addBeardBoxesFromRootsCSV_BBoxFit(const std::string& rootsPath,
                                           double hairLen,
                                           double hairRad,
                                           const Material& hairMat);
    Mesh _addBeardBoxesFromRootsCSV_BBoxFit(const std::string& rootsPath,
                                            double hairLen,
                                            double hairRad,
                                            const Material& hairMat);



    int getTriId() const;

    // Mesh.h にデバッグ用メンバ関数を追加
    void debugTriangleNormals() const {
        if (triangles.empty()) return;

        // メッシュ全体の中心
        Eigen::Vector3d center = Eigen::Vector3d::Zero();
        for (const auto& tri : triangles) {
            center += tri.centroid();
        }
        center /= (double)triangles.size();

        int inwardCount = 0;
        for (int i = 0; i < (int)triangles.size(); ++i) {
            const auto& tri = triangles[i];
            Eigen::Vector3d c = tri.centroid();
            Eigen::Vector3d n = tri.n;

            Eigen::Vector3d outDir = (c - center);       // 中心→三角形の方向
            double dot = outDir.dot(n);

            if (dot < 0.0) {
                // 法線が中心に向かっている（＝内向きっぽい）
                ++inwardCount;
                std::cout << "triangle " << i
                          << " seems inward. dot=" << dot
                          << " normal=" << n.transpose()
                          << " centroid=" << c.transpose()
                          << std::endl;
            }
        }
        std::cout << "Inward-ish triangles: "
                  << inwardCount << " / " << triangles.size() << std::endl;
    }


};



#endif //MESH_H
