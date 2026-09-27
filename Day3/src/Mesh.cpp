//
// Created by nakat on 2025/10/30.
//
#include "Ray.h"
#include "Mesh.h"
#include "Triangle.h"
#include "Image.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

void Mesh::setTriId(int i) const{
    triId = i;
};
int Mesh::getTriId() const {
    return triId;
}

Eigen::Vector3d Mesh::getPos() const {
    return (meshBBox.min+meshBBox.max)/2;
}
Eigen::Vector3d Mesh::getSize() const {
    return meshBBox.max-meshBBox.min;
}

Mesh::Mesh(const std::vector<Triangle> &tri) {

    triangles=tri;
    buildBVH();
}

void Mesh::setTexture(const std::string& fileName) {
    Image img;
    if (!img.loadImage(fileName)) {
        std::cerr << "[Mesh::setTexture] failed to load " << fileName << std::endl;
        return;
    }
    if (img.width <= 0 || img.height <= 0 || triangles.empty()) return;

    auto sampleTex = [&](const Eigen::Vector2d& uv) -> Color {
        double u = uv.x() - std::floor(uv.x());
        double v = uv.y() - std::floor(uv.y());

        double x = u * (img.width  - 1);
        double y = v * (img.height - 1);
        int xx = u * (img.width  - 1);
        int yy = v * (img.height - 1);

        int x0 = std::clamp((int)std::floor(x), 0, img.width  - 1);
        int y0 = std::clamp((int)std::floor(y), 0, img.height - 1);
        int x1 = std::clamp(x0 + 1, 0, img.width  - 1);
        int y1 = std::clamp(y0 + 1, 0, img.height - 1);

        double tx = x - x0;
        double ty = y - y0;

        auto at = [&](int xi, int yi) -> Color {
            return img.pixels[yi * img.width + xi];
        };
        auto lerp = [](const Color& a, const Color& b, double t) {
            return a + (b - a) * t;
        };

        Color c00 = at(x0, y0);
        Color c10 = at(x1, y0);
        Color c01 = at(x0, y1);
        Color c11 = at(x1, y1);

        Color v0 = lerp(c00, c10, tx);
        Color v1 = lerp(c01, c11, tx);
        return lerp(v0, v1, ty);
    };

    for (auto& tri : triangles) {
        const Eigen::Vector2d uv = (tri.uv0+tri.uv1+tri.uv2)/3.0 ;
        const Color texColor = sampleTex(uv);
        tri.material.color = texColor;


    }
}

bool Mesh::loadTexture(const std::string& fileName) {
    if (!m_texture.loadImage(fileName)) {
        std::cerr << "[Mesh::loadTexture] failed to load " << fileName << std::endl;
        m_hasTexture = false;
        return false;
    }
    m_hasTexture = true;
    return true;
}

Color Mesh::getTexture(const Eigen::Vector2d& uv) const {
    if (!m_hasTexture || m_texture.width <= 0 || m_texture.height <= 0) {
        return Color(1.0, 1.0, 1.0);
    }

    double u = uv.x() - std::floor(uv.x());
    double v = uv.y() - std::floor(uv.y());

    double x = u * (m_texture.width  - 1);
    double y = v * (m_texture.height - 1);

    int x0 = std::clamp((int)std::floor(x), 0, m_texture.width  - 1);
    int y0 = std::clamp((int)std::floor(y), 0, m_texture.height - 1);
    int x1 = std::clamp(x0 + 1, 0, m_texture.width  - 1);
    int y1 = std::clamp(y0 + 1, 0, m_texture.height - 1);

    double tx = x - x0;
    double ty = y - y0;

    auto at = [&](int xi, int yi) -> Color {
        return m_texture.pixels[yi * m_texture.width + xi];
    };
    auto lerp = [](const Color& a, const Color& b, double t) {
        return a + (b - a) * t;
    };

    Color c00 = at(x0, y0);
    Color c10 = at(x1, y0);
    Color c01 = at(x0, y1);
    Color c11 = at(x1, y1);

    Color v0 = lerp(c00, c10, tx);
    Color v1 = lerp(c01, c11, tx);
    return lerp(v0, v1, ty);
}
void Mesh::setBSSRDFParams(double scale,double eta,Eigen::Vector3d sigma_a,Eigen::Vector3d sigma_s)
{
    for (Triangle& T:triangles) {
        double kd=T.material.kd;
        Eigen::Vector3d color=T.material.color;



        double emi=T.material.emission;

        T.material=Material(color,kd,emi,eta,sigma_a,sigma_s);
        T.material.scale=scale;

    }


}

void Mesh::buildBVH() {
    // 各三角形のAABBを計算
    // Mesh::buildBVH() の最初など
    const int N = (int)triangles.size();
    triBBox.resize(N);
    primIndices.resize(N);
    for (int i = 0; i < N; ++i) {
        triBBox[i] = triangles[i].computeAABB(); // Triangle側にAABB計算を用意
        primIndices[i] = i;                       // 0,1,2,...,N-1
    }
    bvhRoot = buildRecursive(0, N, 0);


    // BVHを再帰的に構築
    bvhRoot = buildRecursive(0, triangles.size(), 0);




    std::cout << "try size: " << triangles.size() << std::endl;
    std::cout<<"built BVH :BVH size ="<< triangleBBoxes.size()<<std::endl;

    std::cout << "root range=" << bvhRoot->range
          << " isLeaf=" << bvhRoot->isLeaf << "\n";
    std::cout << "L range=" << bvhRoot->left->range
              << "  R range=" << bvhRoot->right->range << "\n";
    std::cout << "root bbox min=" << bvhRoot->bbox.min.transpose()
              << " max=" << bvhRoot->bbox.max.transpose() << "\n"
              << "left bbox min=" << bvhRoot->left->bbox.min.transpose()
              << " max=" << bvhRoot->left->bbox.max.transpose() << "\n"
              << "right bbox min=" << bvhRoot->right->bbox.min.transpose()
              << " max=" << bvhRoot->right->bbox.max.transpose() << "\n";


    buildBoundsOnly();
}

BVHNode* Mesh::buildRecursive(int start, int end, int depth) {
    BVHNode* node = new BVHNode();
    node->start = start;
    node->range = end - start;

    // --- (1) primIndices でノードAABBを作る ---
    AABB box;                    // 無効箱で初期化される実装にしておく
    for (int i = start; i < end; ++i) {
        const int id = primIndices[i];
        box = AABB::merge(box, triBBox[id]);
    }
    node->bbox = box;

    // 葉条件
    const int LeafSize = 4;
    if (node->range <= LeafSize) { node->isLeaf = true; return node; }

    // --- (2) 最大エクステント軸で分割し、primIndices を nth_element ---
    Eigen::Vector3d ext = node->bbox.max - node->bbox.min;
    int axis = (ext[1] > ext[0]) ? 1 : 0;
    if (ext[2] > ext[axis]) axis = 2;

    const int mid = (start + end) / 2;
    std::nth_element(primIndices.begin() + start,
                     primIndices.begin() + mid,
                     primIndices.begin() + end,
                     [&](int ia, int ib){
                         return triangles[ia].centroid()[axis] < triangles[ib].centroid()[axis];
                     });

    node->left  = buildRecursive(start, mid,  depth + 1);
    node->right = buildRecursive(mid,   end,  depth + 1);
    return node;
}

//デバッグ
void Mesh::buildBoundsOnly() {
    triangleBBoxes.resize(triangles.size());
    // 三角形AABBとメッシュAABB
    bool first = true;
    for (int i = 0; i < (int)triangles.size(); ++i) {
        triangleBBoxes[i] = triangles[i].computeAABB();
        if (first) { meshBBox = triangleBBoxes[i]; first = false; }
        else       { meshBBox = AABB::merge(meshBBox, triangleBBoxes[i]); }
    }
}
bool Mesh::hitSingleBBox(const Ray& ray, RayHit& hit) const {
    // 既に見つけた最近ヒット距離を上限に枝刈り
    const double tmin = 1e-4;
    if (!meshBBox.hit(ray, tmin, hit.t)) return false;

    bool any = false;
    for (int i = 0; i < (int)triangles.size(); ++i) {
        RayHit tmp;
        if (triangles[i].hit(ray, tmp) && tmp.t > tmin && tmp.t < hit.t) {
            hit = tmp;
            setTriId(i);
            any = true;
        }
    }
    return any;
}
//kokomade


bool Mesh::hitBVH(const Ray& ray, RayHit& hit, const BVHNode* node) const {

    if (!node || !node->bbox.hit(ray, 1e-4, hit.t)) return false;
    const double tmin = 1e-4;
    bool any = false;
    if (node->isLeaf) {
        for (int i = node->start; i < node->start + node->range; ++i) {
            const int triId = primIndices[i];      // ← 直接 triangles[i] としない
            RayHit tmp;
            if (triangles[triId].hit(ray, tmp) && tmp.t < hit.t) {
                hit = tmp;
                setTriId(triId);
                any = true;
            }
        }
        return any;
    }
    // 近い子から（AABBの入りtで判定）
    auto entryT = [&](const BVHNode* nd)->double{
        // 簡易：レイとAABBのt範囲を再計算し entry を返す
        double t0 = tmin, t1 = hit.t;
        for (int a=0;a<3;++a){
            double invD = 1.0 / ray.dir[a];
            double ta = (nd->bbox.min[a] - ray.org[a])*invD;
            double tb = (nd->bbox.max[a] - ray.org[a])*invD;
            if (invD < 0.0) std::swap(ta,tb);
            t0 = std::max(t0, ta);
            t1 = std::min(t1, tb);
            if (t1 <= t0) return DBL_MAX;
        }
        return t0;
    };

    const double tL = entryT(node->left);
    const double tR = entryT(node->right);
    const BVHNode* first  = (tL < tR ? node->left  : node->right);
    const BVHNode* second = (tL < tR ? node->right : node->left);

    any |= hitBVH(ray, hit, first);
    // 1個目で hit.t が小さくなったら、2個目はより強く枝刈りされる
    any |= hitBVH(ray, hit, second);
    return any;
}
bool Mesh::hitBVH(const Ray& ray, RayHit& hit,
                  const BVHNode* node, int ignoreTriId) const {
    if (!node || !node->bbox.hit(ray, 1e-4, hit.t)) {
    ///    if (!node->bbox.hit(ray, 1e-4, hit.t)) {
    ///        std::cout<<"bbox.min="<<node->bbox.min<<std::endl;
    ///        std::cout<<"bbox.max="<<node->bbox.max<<std::endl;
    ///        std::cout<<"bbox.hit()="<<node->bbox.hit(ray, 1e-4, hit.t)<<std::endl;
    ///    }
        return false;
    }
    const double tmin = 1e-4;
    bool any = false;

    if (node->isLeaf) {
        //std::cout<<"hit"<<std::endl;
        for (int i = node->start; i < node->start + node->range; ++i) {
            const int triId = primIndices[i];

            // ★ 自己交差を無視
            if (triId == ignoreTriId) continue;

            RayHit tmp;
            if (triangles[triId].hit(ray, tmp) && tmp.t < hit.t) {
                hit = tmp;
                setTriId(triId);   // Mesh 内部に triId を記録する用
                any = true;
            }
        }
        return any;
    }

    auto entryT = [&](const BVHNode* nd)->double{
        double t0 = tmin, t1 = hit.t;
        for (int a=0; a<3; ++a) {
            double invD = 1.0 / ray.dir[a];
            double ta = (nd->bbox.min[a] - ray.org[a]) * invD;
            double tb = (nd->bbox.max[a] - ray.org[a]) * invD;
            if (invD < 0.0) std::swap(ta, tb);
            t0 = std::max(t0, ta);
            t1 = std::min(t1, tb);
            if (t1 <= t0) return DBL_MAX;
        }
        return t0;
    };

    const double tL = entryT(node->left);
    const double tR = entryT(node->right);
    const BVHNode* first  = (tL < tR ? node->left  : node->right);
    const BVHNode* second = (tL < tR ? node->right : node->left);

    any |= hitBVH(ray, hit, first,  ignoreTriId);
    any |= hitBVH(ray, hit, second, ignoreTriId);
    return any;
}




// Mesh.cpp
bool Mesh::hit(const Ray &ray, RayHit &hit, int ignoreTriId) const {
    hit.t = DBL_MAX;

    return hitBVH(ray, hit, bvhRoot, ignoreTriId);
}


Eigen::Vector3d Mesh::getKd() const {
    Triangle T=triangles[triId];

    return T.getMaterial().kd * T.getMaterial().color;
}

Eigen::Vector2d Mesh::getUV(const Eigen::Vector3d point) const{
    Triangle T=triangles[triId];
    return T.getUV(point);
}
bool Mesh::isLight() const {
    Triangle T=triangles[triId];
    return T.getMaterial().emission >0.0;
}

Eigen::Vector3d Mesh::getEmission() const{
    Triangle T= triangles[triId];
    return T.getMaterial().emission * T.getMaterial().color;
}
Material Mesh::getMaterial() const{
    Triangle T=triangles[triId];

    return  T.getMaterial();
}

//VBh

void Mesh::scale(double s) {
    std::vector<Triangle>tmp=triangles;
    for (int i=0; i<tmp.size(); ++i) {
        tmp[i].v0*=s;
        tmp[i].v1*=s;
        tmp[i].v2*=s;
    }
    triangles=tmp;

    // 三角形をスケーリング

    primIndices.clear();
    triBBox.clear();
    bvhRoot=nullptr;

    buildBVH();


}
static inline std::vector<Eigen::Vector2d> loadRootsCSV(const std::string& path) {
    std::ifstream ifs(path);
    std::string line;
    std::vector<Eigen::Vector2d> roots;
    if (!ifs) return roots;

    std::getline(ifs, line); // header
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string sx, sz;
        if (!std::getline(ss, sx, ',')) continue;
        if (!std::getline(ss, sz, ',')) continue;
        roots.emplace_back(std::stod(sx), std::stod(sz));
    }
    return roots;
}

static inline void appendHairBoxTris(
    std::vector<Triangle>& tris,
    const Eigen::Vector3d& center,
    const Eigen::Vector3d& size,      // (sx, sy, sz)
    const Material& mat
){
    Eigen::Vector3d h = 0.5 * size; // half extents

    // 8 corners
    Eigen::Vector3d v000 = center + Eigen::Vector3d(-h.x(), -h.y(), -h.z());
    Eigen::Vector3d v001 = center + Eigen::Vector3d(-h.x(), -h.y(), +h.z());
    Eigen::Vector3d v010 = center + Eigen::Vector3d(+h.x(), -h.y(), -h.z());
    Eigen::Vector3d v011 = center + Eigen::Vector3d(+h.x(), -h.y(), +h.z());
    Eigen::Vector3d v100 = center + Eigen::Vector3d(-h.x(), +h.y(), -h.z());
    Eigen::Vector3d v101 = center + Eigen::Vector3d(-h.x(), +h.y(), +h.z());
    Eigen::Vector3d v110 = center + Eigen::Vector3d(+h.x(), +h.y(), -h.z());
    Eigen::Vector3d v111 = center + Eigen::Vector3d(+h.x(), +h.y(), +h.z());

    auto add = [&](const Eigen::Vector3d& a,
                   const Eigen::Vector3d& b,
                   const Eigen::Vector3d& c){
        Triangle t(a,b,c, mat);   // ← Triangleのコンストラクタに合わせて調整
        tris.push_back(t);
    };

    // 6 faces * 2 tris
    // bottom (-y)
    add(v000, v010, v011); add(v000, v011, v001);
    // top (+y)
    add(v100, v101, v111); add(v100, v111, v110);
    // -x
    add(v000, v001, v101); add(v000, v101, v100);
    // +x
    add(v010, v110, v111); add(v010, v111, v011);
    // -z
    add(v000, v100, v110); add(v000, v110, v010);
    // +z
    add(v001, v011, v111); add(v001, v111, v101);
}


void Mesh::addBeardFromRootsCSV(const std::string& rootsPath,
                               double hairLen,
                               double hairRad,
                               double ySurface,
                               bool rootsAreMeters,
                               const Material& hairMat)
{
    auto roots = loadRootsCSV(rootsPath);

    const double unit = rootsAreMeters ? 1000.0 : 0.1; // 自作がmmなら meter->mm
    for (auto& rz : roots) {
        double x = rz.x() * unit;
        double z = rz.y() * unit;

        Eigen::Vector3d size(2.0*hairRad, hairLen, 2.0*hairRad);
        Eigen::Vector3d center(x, ySurface + hairLen*0.5, z);

        appendHairBoxTris(this->triangles, center, size, hairMat);
    }

    // trianglesを増やしたのでBVHを作り直す
    this->buildBVH();
}
struct RootsData {
    std::vector<Eigen::Vector2d> roots; // (x,z)
    double minX=1e30, maxX=-1e30;
    double minZ=1e30, maxZ=-1e30;
};

static RootsData loadRootsCSV_withMinMax(const std::string& path) {
    RootsData rd;
    std::ifstream ifs(path);
    std::string line;
    if (!ifs) return rd;

    std::getline(ifs, line); // header
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string sx, sz;
        if (!std::getline(ss, sx, ',')) continue;
        if (!std::getline(ss, sz, ',')) continue;

        double x = std::stod(sx);
        double z = std::stod(sz);

        rd.roots.emplace_back(x, z);
        rd.minX = std::min(rd.minX, x); rd.maxX = std::max(rd.maxX, x);
        rd.minZ = std::min(rd.minZ, z); rd.maxZ = std::max(rd.maxZ, z);
    }
    return rd;
}

static inline double mapLinear(double v, double a0, double a1, double b0, double b1) {
    if (std::abs(a1 - a0) < 1e-12) return 0.5 * (b0 + b1); // 全部同値の保険
    double t = (v - a0) / (a1 - a0);
    return b0 + t * (b1 - b0);
}
void Mesh::addBeardBoxesFromRootsCSV_BBoxFit(const std::string& rootsPath,
                                             double hairLen,
                                             double hairRad,
                                             const Material& hairMat)
{

    RootsData rd = loadRootsCSV_withMinMax(rootsPath);
    if (rd.roots.empty()) return;

    // ★ bbox を基準にする（BVH rootがあるならそれを使う）
    AABB bb = this->bvhRoot->bbox;  // ←あなたの実装に合わせて取得
    const double bx0 = bb.min.x(), bx1 = bb.max.x();
    const double bz0 = bb.min.z(), bz1 = bb.max.z();

    // 生やす高さ（平面の上面）
    const double ySurf = bb.max.y();

    for (auto& r : rd.roots) {
        const double x = mapLinear(r.x(), rd.minX, rd.maxX, bx0, bx1);
        const double z = mapLinear(r.y(), rd.minZ, rd.maxZ, bz0, bz1);

        Eigen::Vector3d size(2.0*hairRad, hairLen, 2.0*hairRad);
        Eigen::Vector3d center(x, ySurf + 0.5*hairLen, z);

        appendHairBoxTris(this->triangles, center, size, hairMat);
    }

    this->buildBVH();
}
Mesh Mesh::_addBeardBoxesFromRootsCSV_BBoxFit(const std::string& rootsPath,
                                             double hairLen,
                                             double hairRad,
                                             const Material& hairMat
                                             )
{
    Mesh B=Mesh();
    RootsData rd = loadRootsCSV_withMinMax(rootsPath);
    if (rd.roots.empty()) return B;

    // ★ bbox を基準にする（BVH rootがあるならそれを使う）
    AABB bb = this->bvhRoot->bbox;  // ←あなたの実装に合わせて取得
    const double bx0 = bb.min.x(), bx1 = bb.max.x();
    const double bz0 = bb.min.z(), bz1 = bb.max.z();

    // 生やす高さ（平面の上面）
    const double ySurf = bb.max.y();

    for (auto& r : rd.roots) {
        const double x = mapLinear(r.x(), rd.minX, rd.maxX, bx0, bx1);
        const double z = mapLinear(r.y(), rd.minZ, rd.maxZ, bz0, bz1);

        Eigen::Vector3d size(2.0*hairRad, hairLen, 2.0*hairRad);
        Eigen::Vector3d center(x, ySurf + 0.5*hairLen, z);

        appendHairBoxTris(B.triangles, center, size, hairMat);
    }

    B.buildBVH();

    return B;
}





