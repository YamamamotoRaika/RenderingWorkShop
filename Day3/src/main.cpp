#include <filesystem>
#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <map>
#include "Body.h"
#include "Camera.h"
#include "Renderer.h"
#include "Triangle.h"

struct InternalMaterial
{
    std::string material_name;
    std::string texture_name;
    Eigen::Vector3d kd;
    bool valid;
};

void intersectTest() {
    const Sphere sphere(1, Eigen::Vector3d::Zero());
    const Triangle triangle(Eigen::Vector3d(0,1,0),Eigen::Vector3d(-1,-1,0),Eigen::Vector3d(1,-1,0));
    Ray ray(Eigen::Vector3d(0, 0, 10), Eigen::Vector3d(0, 0, -1));
    RayHit hit;
    sphere.hit(ray, hit);
    std::cout << "sphere" << std::endl;
    std::cout << "t:\t" << hit.t << std::endl;
    std::cout << "normal:\t(" << hit.normal.transpose() << ")" << std::endl;

    triangle.hit(ray, hit);
    std::cout << "triangle" << std::endl;
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
    const unsigned int samples = 10000;
    const auto image1 = renderer.render();
    //const auto image2 = renderer.directIlluminationRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();

    image1.save("sample_image1.png");
    //image2.apply_reinhard_extended_tone_mapping().save("sample.png");
}
void Bokasi(const std::string fileName) {

    Image img;  // サイズは loadImage の中で決まる

    if (!img.loadImage(fileName)) {
        std::cerr << "[Renderer::loadTex] failed to load " << fileName << std::endl;
        return;
    }

   img.apply_reinhard_extended_tone_mapping().apply_gaussian_blur(7,7,2.0);

    img.save("sss_blurred.png");

}

bool loadMtl(std::map<std::string,Material> &materials ,const std::string &filename) {
    //file読み込みのための準備
    FILE* f= NULL;

    const int BUFFER_SIZE = 4096;
    char line[BUFFER_SIZE];
    char material_name[BUFFER_SIZE];
    f = fopen(filename.c_str(), "r");

    if (!f) {
        std::cout << "MTL File not found" << std::endl;
        //std::filesystem::path cwd = std::filesystem::current_path();
        //std::cout << "現在の作業ディレクトリ: " << cwd << std::endl;
        return false;
    }

    Material mat=Material();
    mat.color << 1.0, 1.0, 1.0;
    mat.kd=1;
    bool valid = false;

    while (fgets(line, BUFFER_SIZE, f)!=NULL) {
        if ( strncmp( line, "newmtl", 6 ) == 0 )
        {
            if( valid ) {
                materials[material_name] = mat;
                valid = false;
            }

            sscanf( line, "newmtl %s", material_name );
            if ( strlen( material_name ) == 0 )
                continue;


        //mtlが切り替わったので値を初期化している
            mat.color << 0.55, 0.38, 0.32;
            mat.kd=1;
            valid = true;
        }
        // 先頭の文字を識別
        if ( line[0] == 'K' )
        {
            //  (Kd)はColorにぶち込む
            if ( line[1] == 'd' )
            {
                float r, g, b;
                sscanf( line, "Kd %f %f %f", &r, &g, &b );

                mat.color << 0.55, 0.38, 0.32;
                mat.kd = (r+g+b)/3;
            }
        }

        if ( strncmp( line, "map_Kd", 6 ) == 0 )
        {/*
            sscanf( line, "map_Kd %s", texture_name );

            if ( strlen( texture_name ) == 0 )
                continue;

            material.texture_name = texture_name;
            */
        }

    }

    if( valid )
    {
        materials[material_name]= mat ;
        valid = false;
    }

    fclose( f );
    return true;
}







bool loadObj(  std::vector<Body> &model , const std::string &filename ) {
    //file読み込みのための準備
    FILE* f= NULL;

    const int BUFFER_SIZE = 4096;
    char line[BUFFER_SIZE];
    char material_name[BUFFER_SIZE];
    f = fopen(filename.c_str(), "r");

    if (!f) {
        std::cout << "OBJ File not found" << std::endl;
        //std::filesystem::path cwd = std::filesystem::current_path();
        //std::cout << "現在の作業ディレクトリ: " << cwd << std::endl;
        return false;
    }
    //必要な要素を宣言、いったんmtlファイルの読み込みは考えないで実装,


    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vertices;
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vertex_normals;
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> tex_coords;

    Eigen::Vector3d kd;
    Material mat=Material(codeToColor("#e597b2"), 1.0);
    std::map<std::string,Material> materials ;


    //objfileの中身を一行ずつ読み込んでいく
    while( fgets( line, BUFFER_SIZE, f ) != NULL ) {


        //mtllibの時、マテリアル のファイルを読み込む
        if ( strncmp( line, "mtllib", 6 ) == 0 )
        {
            sscanf( line, "mtllib %s", material_name );

            if ( strlen( material_name ) > 0 )
            {
                if( !loadMtl( materials,material_name  ) )
                {
                    std::cout << "Could not read mtl file: " << material_name << std::endl;
                    fclose( f );
                    return false;
                }

            }
        }

        // マテリアルを使う。いったん中身はコメントアウト
        if ( strncmp( line, "usemtl", 6 ) == 0 )
        {
            sscanf( line, "usemtl %s", material_name );

            /*
            if( triangles.triangles.size() > 0 )
            {
                out_internal_triangles.push_back( triangles );
                resetInternalTriangles( triangles );
            }
            */
            std::cout << material_name << std::endl;
            mat=materials[material_name];
        }

        //頂点情報の処理
        if (line[0]=='v') {
            //法線情報の処理
            if (line[1]=='n') {
                float x, y, z;
                sscanf( line, "vn %f %f %f", &x, &y, &z );
                vertex_normals.emplace_back(x,y,z);
            }
            //テクスチャ座標の処理
            else if (line[1]=='t') {
                float u, v;
                sscanf( line, "vt %f %f", &u, &v );
                tex_coords.emplace_back( u, 1.0 - v );
            }
            else {
                float x, y, z;
                sscanf( line, "v %f %f %f", &x, &y, &z );
                vertices.emplace_back(x,y,z);
            }



        }

        //ポリゴン情報（面）の処理
        if (line[0]=='f') {
            char* tp = strtok( &(line[2]), " " );
            std::vector<Eigen::Vector3i> vertices_number;

            std::cout<<mat.kd<<std::endl;
            while( tp != NULL )
            {
                int vid = -1, vtid = -1, vnid = -1;
                sscanf( tp, "%d/%d/%d", &vid, &vtid, &vnid );

                // objファイルの中身は1から始まるので、１引く
                vertices_number.emplace_back( vid-1, vtid-1, vnid-1 );

                tp = strtok( NULL, " " );
            }

            //取得した頂点情報を基に三角形を作成する
            for (int i=0;i<vertices_number.size()-2;i++) {
                //法線情報も考えないよ,いったんね
                Triangle t=Triangle(vertices[vertices_number[0][0]],
                                    vertices[vertices_number[i+1][0]],
                                    vertices[vertices_number[i+2][0]]);

                Body body(t, mat);
                model.push_back(body);
            }
        }

        if (model.size()>0) {

        }





    }


    fclose(f);
    std::cout << "Loaded Model: " << filename << std::endl;
    return true;
}



bool _loadObj(std::vector<Triangle> &model, const std::string &filename) {
    FILE* f = fopen(filename.c_str(), "r");
    if (!f) {
        std::cout << "File not found: " << filename << std::endl;
        return false;
    }

    const int BUFFER_SIZE = 4096;
    char line[BUFFER_SIZE];
    char material_name[BUFFER_SIZE];

    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vertices;
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vertex_normals;
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> tex_coords;

    Material mat = Material(Color(  0.55, 0.38, 0.32), 1.0);
    std::map<std::string, Material> materials;

    // OBJのあるディレクトリ
    std::filesystem::path objPath(filename);
    std::filesystem::path baseDir = objPath.parent_path();

    auto safeUV = [&](int uvIdx)->Eigen::Vector2d {
        if (uvIdx >= 0 && uvIdx < (int)tex_coords.size()) return tex_coords[uvIdx];
        return Eigen::Vector2d(0.0, 0.0); // vtが無いOBJ用
    };

    while (fgets(line, BUFFER_SIZE, f) != NULL) {
        // 行頭の空白・タブを飛ばす
        char* p = line;
        while (*p == ' ' || *p == '\t') ++p;

        // コメント/空行
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') continue;

        // mtllib
        if (strncmp(p, "mtllib", 6) == 0) {
            if (sscanf(p, "mtllib %s", material_name) == 1) {
                std::filesystem::path mtlPath = baseDir / material_name; // ./apple.mtl もOK
                if (!loadMtl(materials, mtlPath.string())) {
                    std::cout << "[Warn] Could not read mtl file: " << mtlPath.string()
                              << " (continue with default material)\n";
                    // ここでreturnしない
                }
            }
            continue;
        }

        // usemtl
        if (strncmp(p, "usemtl", 6) == 0) {
            if (sscanf(p, "usemtl %s", material_name) == 1) {
                auto it = materials.find(material_name);
                if (it != materials.end()) {
                    mat = it->second;
                } else {
                    std::cout << "[Warn] usemtl not found in materials: " << material_name
                              << " (keep current material)\n";
                }
            }
            continue;
        }

        // v / vt / vn
        if (p[0] == 'v') {
            if (p[1] == 'n') {
                float x, y, z;
                if (sscanf(p, "vn %f %f %f", &x, &y, &z) == 3)
                    vertex_normals.emplace_back(x, y, z);
            } else if (p[1] == 't') {
                float u, v;
                if (sscanf(p, "vt %f %f", &u, &v) == 2)
                    tex_coords.emplace_back(u, 1.0 - v);
            } else {
                float x, y, z;
                if (sscanf(p, "v %f %f %f", &x, &y, &z) == 3)
                    vertices.emplace_back(x, y, z);
            }
            continue;
        }

        // f
        if (p[0] == 'f') {
            // 区切りに \t\r\n も入れる（最後のトークンに改行が残る事故を防ぐ）
            char* tp = strtok(p + 1, " \t\r\n");
            std::vector<Eigen::Vector3i> idx; // (vid, vtid, vnid)

            while (tp != NULL) {
                int vid = -1, vtid = -1, vnid = -1;

                if (strstr(tp, "//")) {
                    sscanf(tp, "%d//%d", &vid, &vnid); // v//vn
                } else {
                    char* s1 = strchr(tp, '/');
                    char* s2 = s1 ? strchr(s1 + 1, '/') : nullptr;

                    if (!s1)            sscanf(tp, "%d", &vid);              // v
                    else if (!s2)       sscanf(tp, "%d/%d", &vid, &vtid);     // v/vt
                    else                sscanf(tp, "%d/%d/%d", &vid, &vtid, &vnid); // v/vt/vn
                }

                idx.emplace_back(vid - 1, vtid - 1, vnid - 1);
                tp = strtok(NULL, " \t\r\n");
            }

            // fan triangulation
            for (int i = 0; i + 2 < (int)idx.size(); ++i) {
                const auto &i0 = idx[0];
                const auto &i1 = idx[i + 1];
                const auto &i2 = idx[i + 2];

                // 頂点Indexの安全チェック（最低限）
                if (i0[0] < 0 || i1[0] < 0 || i2[0] < 0 ||
                    i0[0] >= (int)vertices.size() ||
                    i1[0] >= (int)vertices.size() ||
                    i2[0] >= (int)vertices.size()) {
                    continue;
                }

                const Eigen::Vector3d &p0 = vertices[i0[0]];
                const Eigen::Vector3d &p1 = vertices[i1[0]];
                const Eigen::Vector3d &p2 = vertices[i2[0]];

                Eigen::Vector2d u0 = safeUV(i0[1]);
                Eigen::Vector2d u1 = safeUV(i1[1]);
                Eigen::Vector2d u2 = safeUV(i2[1]);

                Eigen::Vector3d faceN = (p1 - p0).cross(p2 - p0);
                if (faceN.squaredNorm() == 0.0) continue;
                faceN.normalize();

                Triangle t(p0, p1, p2, u0, u1, u2, faceN);
                t.material = mat;
                model.push_back(t);
            }
            continue;
        }
    }

    fclose(f);
    std::cout << "Loaded Model: " << filename << " (triangles=" << model.size() << ")\n";
    return true;
}

//Obj型の導入に際して
void ObjTest() {



    std::vector<Triangle> Meshes {
    };

    const std::vector<Body> lights {
       //Body(Sphere(1, Eigen::Vector3d(15, 15, 15)), Material(codeToColor("#ffffff"), 1.0, 300)),
        Body(Sphere(1, Eigen::Vector3d(0, 10, 0)), Material(codeToColor("#ffffff"), 1.0, 10)),
        Body(Sphere(1, Eigen::Vector3d(-20, 0, -10)), Material(codeToColor("#ffffff"), 1.0, 10)),
        Body(Sphere(1, Eigen::Vector3d(-1.7, -1.2, 1.1)), Material(codeToColor("#ffffff"), 1.0, 1)),


        //Body(Sphere(5, Eigen::Vector3d(2, 10, 10)), Material(codeToColor("#ffffff"), 1.0, 300))
};

    _loadObj(Meshes,"../Day3/lpshead/head.OBJ");
    //変更点　まずMeshsesをそのまま入れる→
    //Meshクラスにしたら？→render側でマテリアルを直接参照している部分を関数呼び出しにしたら行けた。
    Mesh M=Mesh(Meshes);
    Mesh B=Mesh();


    double room_r=50;

    Eigen::Vector3d S = M.getSize();


    //double size=20/S[1];
    M.scale(5);;
    M.setBSSRDFParams(1,1.3,Eigen::Vector3d(0.0011, 0.0024, 0.014),Eigen::Vector3d(0.74, 0.88, 1.01));
    //texture適用　
    M.setTexture("./lpshead/lambertian.jpg");




    std::vector<Body> bodies {

        M
    };
    //std::cout << "bodies.size() = " << bodies.size() << std::endl;
    std::cout << " = " << M.getSize() << std::endl;





        for(const auto & light : lights) {
            bodies.push_back(light);
        }

    Eigen::Vector3d x = M.getPos();


    const Eigen::Vector3d campos(-1.5, -1.0, 1.3);  // カメラを少し斜め上から
    const Eigen::Vector3d lookat(0.5,-0.2,0); // モデル中心に向ける

    //const Eigen::Vector3d campos(0, 0, 2.3);  // カメラを少し斜め上から
    //const Eigen::Vector3d lookat(0,-0.2,0); // モデル中心に向ける
    const Eigen::Vector3d camdir = lookat - campos;
    const double fov = 45;                    // 垂直FoV
    const double aspect = 4.0 / 3.0;
    const int height = 540;

    const Camera camera(campos, camdir, height, aspect, fov);




    /// 背景色はわかりやすく灰色
    const Renderer renderer(bodies, camera, Color(0.1, 0.1, 0.1));
   //const auto image = renderer.render().apply_reinhard_extended_tone_mapping().apply_gamma_correction();
    //image.save("TEST.png");
    const unsigned int samples = 1000;
    //const auto image = renderer.directIlluminationRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();
    //image.save("BRDF.png");


    //const auto image = renderer._SSSdirectIlluminationRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();
    //image.save("BSSRDF_Dipole.png");

    const auto image = renderer.KAI_SSSdirectIlluminationRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();

    //image.save("gazo_applied.png");
    //image.save("gazo_nobeard.png");
    //image.save("gazo_fullbeard.png");

    image.save("face_notexture_applied.png");
    //image.save("face_notexture_nobeard.png");
    //image.save("face_notexture_fullbeard.png");

    //image.save("face_texture_applied.png");
    //image.save("face_texture_nobeard_test.png");

    //image.save("face_texture_fullbeard.png");

    // Rd(r) profile (beard-mixed) at camera-center hit
    const bool measureRd = false;
    if (measureRd) {
        const size_t angularSamples = 256;
        const double rMax = 2;
        const int bins = 200;
        renderer.ProfileRadialBSSRDF_BeardMixed_CameraCenter(
            angularSamples, rMax, bins, "Rd_profile.csv");
    }



    //image.save("BSSRDF_Multipole.png");



}
void SSSTest() {
const auto room_r = 120;
    const auto floor_color = codeToColor("#fedcbd");
    //座標修正
    const std::vector<Body> room_walls {
    ///        Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(1,0,0)),
    ///            Material(codeToColor("#2f5d50"), 0.8, 0.0,1.3,Eigen::Vector3d(0.03,0.17,0.48),Eigen::Vector3d(0.74,0.88,1.01),true
    ///            )),
    ///        Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(-1,0,0)),
    ///        Material(codeToColor("#00a3af"), 0.8, 0.0,1.3,Eigen::Vector3d(0.03,0.17,0.48),Eigen::Vector3d(0.74,0.88,1.01),true
    ///    )),
            Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(0,0,-1)),
            Material(floor_color, 0.8, 0.0,1.3,Eigen::Vector3d(0.0011, 0.0024, 0.014),Eigen::Vector3d(0.74, 0.88, 1.01),true
                )),
    ///        Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(0,0.9,0)),
    ///        Material(floor_color, 0.8, 0.0,1.3,Eigen::Vector3d(0.03,0.17,0.48),Eigen::Vector3d(0.74,0.88,1.01),true
    ///           )),
    ///       Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(0,-0.8,0)),
    ///       Material(floor_color, 0.8, 0.0,1.3,Eigen::Vector3d(0.03,0.17,0.48),Eigen::Vector3d(0.74,0.88,1.01),true
    ///           )),

};

std::vector<Body> bodies{
    ///
    ///Body(Box(16.5* Eigen::Vector3d(1,1,1), Eigen::Vector3d(25, -14.5, -10)), Material(Color(0.75, 0.25, 0.25), 0.8, 0.0)),
    ///Body(Box(16.5* Eigen::Vector3d(1,1,1), Eigen::Vector3d(-23, -14.5, 7)), Material(Color(0.99, 0.99, 0.99), 0.8, 0.0)),
    };

    const std::vector<Body> lights {
            Body(Sphere(50, Eigen::Vector3d(0, 34.8, 0)), Material(codeToColor("#e597b2"), 1.0, 30))
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

    const unsigned int samples = 100;
    const auto image = renderer.SSSdirectIlluminationRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();
   // const auto image = renderer._directIlluminationRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();
    //const auto image = renderer.ReferenceSSSRandomWalkRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();

    //image.save("Dipole.png");
    //image.save("BRDF.png");
    //image.save("multipole__200m.png");
    image.save("kensyouyou.png");

 //renderer.Kensyou(samples);


}


int main() {
    std::cout << "Hello, World!" << std::endl;

    intersectTest();
    //実行時間計測のため
    auto start = std::chrono::system_clock::now();


    //Bokasi("../Day3/Test_kougou_5.png");
    ObjTest();
    //SSSTest();
    //roomRenderingSample_Box();
    auto end = std::chrono::system_clock::now();

    // end - start をミリ秒単位で計算する
    std::chrono::duration<double, std::milli> elapsed = end - start;

    // end - start を秒単位で計算する
    std::chrono::duration<double> elapsed2 = end - start;

    // 結果をコンソールに出力する
    std::cout << elapsed.count() << "ms" << std::endl;
    std::cout << elapsed2.count() << "s" << std::endl;

    return 0;
}
