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
//三角形テスト用
void TriangleTest() {
    /// bodiesに光源を追加
    const std::vector<Body> bodies =  {
        Body(Triangle (Eigen::Vector3d(0,1,0),Eigen::Vector3d(-1,-1,0),Eigen::Vector3d(1,-1,0)),Material(Color(0.1, 1, 0.1), 0.8)),
        Body(Triangle (Eigen::Vector3d(3,1,0),Eigen::Vector3d(2,-1,-1),Eigen::Vector3d(4,-1,-1)),Material(Color(0.1, 1, 0.1), 0.8)),
        Body(Triangle (Eigen::Vector3d(-3,1,0),Eigen::Vector3d(-4,-1,1),Eigen::Vector3d(-2,-1,1)),Material(Color(0.1, 1, 0.1), 0.8)),

        Body(Sphere(2.0, Eigen::Vector3d(0, 10, 10)), Material(Color(1, 1, 1), 0.8, 20)),
};

    const Eigen::Vector3d campos(0, 10, 100);
    const Eigen::Vector3d camdir = Eigen::Vector3d(0, 0, 0) - campos;

    const Camera camera(campos, camdir, 320, 16.0 / 9.0, 5);

    /// 背景色はわかりやすく灰色
    const Renderer renderer(bodies, camera, Color(0.1, 0.1, 0.1));
    const unsigned int samples = 5000;
    const auto image1 = renderer.render();
    const auto image2 = renderer.directIlluminationRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();

    image1.save("Triangle_1.png");
    image2.apply_reinhard_extended_tone_mapping().save("Triangle_2.png");
}
//直方体テスト用
void BoxTest() {
    /// bodiesに光源を追加
    const std::vector<Body> bodies =  {
        Body(Box (Eigen::Vector3d(0,0,0),Eigen::Vector3d(1,1,1)),Material(Color(0.1, 1, 0.1), 0.8)),

        Body(Sphere(2.0, Eigen::Vector3d(0, 10, 10)), Material(Color(1, 1, 1), 0.8, 10)),
};

    const Eigen::Vector3d campos(0, 10, -100);
    const Eigen::Vector3d camdir = Eigen::Vector3d(0, 0, 0) - campos;

    const Camera camera(campos, camdir, 320, 9.0 / 16.0, 5);

    /// 背景色はわかりやすく灰色
    const Renderer renderer(bodies, camera, Color(0.1, 0.1, 0.1));
    const unsigned int samples = 5000;
    //const auto image1 = renderer.render();
    const auto image2 = renderer.directIlluminationRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();

    //image1.save("Box_1.png");
     image2.apply_reinhard_extended_tone_mapping().save("Box_2.png");
}

void roomRenderingSample() {
    const auto room_r = 1e5;
    const auto floor_color = codeToColor("#f9c89b");
    const std::vector<Body> room_walls {
            Body(Sphere(room_r, (room_r - 50) * Eigen::Vector3d::UnitX()), Material(codeToColor("#2f5d50"), 0.8, 0.0)),
            Body(Sphere(room_r, -(room_r - 50) * Eigen::Vector3d::UnitX()), Material(codeToColor("#00a3af"), 0.8, 0.0)),
            Body(Sphere(room_r, (room_r - 30) * Eigen::Vector3d::UnitY()), Material(floor_color, 0.8, 0.0)),
            Body(Sphere(room_r, -(room_r - 40) * Eigen::Vector3d::UnitY()), Material(floor_color, 0.8, 0.0)),
            Body(Sphere(room_r, (room_r - 30) * Eigen::Vector3d::UnitZ()), Material(floor_color, 0.8, 0.0)),
    };

    std::vector<Body> bodies {
            ///Body(Sphere(16.5, Eigen::Vector3d(25, -14.5, -10)), Material(Color(0.75, 0.25, 0.25), 0.8, 0.0)),
            //Body(Sphere(16.5, Eigen::Vector3d(-23, -14.5, 7)), Material(Color(0.99, 0.99, 0.99), 0.8, 0.0)),
    };

    const std::vector<Body> lights {
            Body(Sphere(5, Eigen::Vector3d(0, 34.8, 0)), Material(codeToColor("#e597b2"), 1.0, 30))
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
    //const auto image = renderer.render().apply_reinhard_extended_tone_mapping().apply_gamma_correction();

    const unsigned int samples = 1e3;

    //image.save("sample_image.png");

}
//Box型の導入に際して
void roomRenderingSample_Box() {
    const auto room_r = 120;
    const auto floor_color = codeToColor("#f9c89b");
    //座標修正
    const std::vector<Body> room_walls {
      //
      //     Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(1,0,0)), Material(codeToColor("#2f5d50"), 0.8, 0.0)),
      //      Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(-1,0,0)), Material(codeToColor("#00a3af"), 0.8, 0.0)),
            Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(0,0,-1)), Material(floor_color, 0.8, 0.0)),
      //      Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(0,0.9,0)), Material(floor_color, 0.8, 0.0)),
      //      Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(0,-0.8,0)), Material(floor_color, 0.8, 0.0)),

    };

    std::vector<Body> bodies {
       //     Body(Box(16.5* Eigen::Vector3d(1,1,1), Eigen::Vector3d(25, -14.5, -10)), Material(Color(0.75, 0.25, 0.25), 0.8, 0.0)),
       //     Body(Box(16.5* Eigen::Vector3d(1,1,1), Eigen::Vector3d(-23, -14.5, 7)), Material(Color(0.99, 0.99, 0.99), 0.8, 0.0)),
    };

    const std::vector<Body> lights {
            Body(Sphere(5, Eigen::Vector3d(0, 34.8, 0)), Material(codeToColor("#e597b2"), 1.0, 30))
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
    //const auto image = renderer.render().apply_reinhard_extended_tone_mapping().apply_gamma_correction();

    const unsigned int samples = 1e4;
    const auto image = renderer.directIlluminationRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();


    image.save("BoxRoom.png");

}
bool loadMtl(std::map<std::string,Material> &materials ,const std::string &filename) {
    //file読み込みのための準備
    FILE* f= NULL;

    const int BUFFER_SIZE = 4096;
    char line[BUFFER_SIZE];
    char material_name[BUFFER_SIZE];
    f = fopen(filename.c_str(), "r");

    if (!f) {
        std::cout << "File not found" << std::endl;
        //std::filesystem::path cwd = std::filesystem::current_path();
        //std::cout << "現在の作業ディレクトリ: " << cwd << std::endl;
        return false;
    }

    Material mat=Material();
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
            mat.color << 1.0, 1.0, 1.0;
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

                mat.color << r, g, b;
                mat.kd = 1;
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
        std::cout << "File not found" << std::endl;
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


//Obj型の導入に際して
void ObjTest() {


    std::vector<Body> bodies {
    };
/*
    const std::vector<Body> lights {
        Body(Sphere(5, Eigen::Vector3d(0, 34.8, 0)), Material(codeToColor("#e597b2"), 1.0, 30))
};
*/
    loadObj(bodies,"../Day2/apple.obj");


/*
        for(const auto & light : lights) {
            bodies.push_back(light);
        }
*/

    const Eigen::Vector3d campos(0, 40, 160);  // カメラを少し斜め上から
    const Eigen::Vector3d lookat(0, 40, -33); // モデル中心に向ける
    const Eigen::Vector3d camdir = lookat - campos;
    const double fov = 45;                    // 垂直FoV
    const double aspect = 4.0 / 3.0;
    const int height = 540;

    const Camera camera(campos, camdir, height, aspect, fov);




    /// 背景色はわかりやすく灰色
    const Renderer renderer(bodies, camera, Color(0.1, 0.1, 0.1));
    const auto image = renderer.render().apply_reinhard_extended_tone_mapping().apply_gamma_correction();

    const unsigned int samples = 500;
    //const auto image = renderer.directIlluminationRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();


    image.save("ObjTest.png");

}
void SSSTest() {
const auto room_r = 120;
    const auto floor_color = codeToColor("#f9c89b");
    //座標修正
    const std::vector<Body> room_walls {
            //Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(1,0,0)),
             //   Material(codeToColor("#2f5d50"), 0.8, 0.0,1.3,Eigen::Vector3d(0.03,0.17,0.48),Eigen::Vector3d(0.74,0.88,1.01),true
              //  )),
            ///Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(-1,0,0)),
            ///Material(codeToColor("#00a3af"), 0.8, 0.0,1.3,Eigen::Vector3d(0.03,0.17,0.48),Eigen::Vector3d(0.74,0.88,1.01),true
        ///)),
            Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(0,0,-1)),
            Material(floor_color, 0.8, 0.0,1.3,Eigen::Vector3d(0.001, 0.002, 0.003),Eigen::Vector3d(1.0, 1.0, 1.0),true
                )),
          ///  Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(0,0.9,0)),
          ///  Material(floor_color, 0.8, 0.0,1.3,Eigen::Vector3d(0.03,0.17,0.48),Eigen::Vector3d(0.74,0.88,1.01),true
           ///     )),
           /// Body(Box((room_r)* Eigen::Vector3d(1,1,1), room_r * Eigen::Vector3d(0,-0.8,0)),
           /// Material(floor_color, 0.8, 0.0,1.3,Eigen::Vector3d(0.03,0.17,0.48),Eigen::Vector3d(0.74,0.88,1.01),true
           ///     )),

};

std::vector<Body> bodies{
    ///
    ///Body(Box(16.5* Eigen::Vector3d(1,1,1), Eigen::Vector3d(25, -14.5, -10)), Material(Color(0.75, 0.25, 0.25), 0.8, 0.0)),
    ///Body(Box(16.5* Eigen::Vector3d(1,1,1), Eigen::Vector3d(-23, -14.5, 7)), Material(Color(0.99, 0.99, 0.99), 0.8, 0.0)),
    };

    const std::vector<Body> lights {
            Body(Sphere(5, Eigen::Vector3d(0, 34.8, 0)), Material(codeToColor("#e597b2"), 1.0, 30))
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

    const unsigned int samples = 500;
    const auto image = renderer.SSSdirectIlluminationRender(samples).apply_reinhard_extended_tone_mapping().apply_gamma_correction();

    image.save("SSSTest.png");



}

int main() {
    std::cout << "Hello, World!" << std::endl;

    intersectTest();
    //実行時間計測のため
    auto start = std::chrono::system_clock::now();


    SSSTest();
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
