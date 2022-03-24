#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include <algorithm>
const int width  = 800;
const int height = 800;
const int depth  = 255;

Model *model = NULL;
int *zbuffer = NULL;
Vec3f light_dir = Vec3f(0,-1,-1).normalize();
//摄像机位置
Vec3f eye(2,1,3);
//焦点位置
Vec3f center(0,0,1);

//视角矩阵，用于将(-1,1),(-1,1),(-1,1)映射到(1/8w,7/8w),(1/8h,7/8h),(0,255)
Matrix viewport(int x, int y, int w, int h) {
    Matrix m = Matrix::identity(4);
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
    return m;
}

//朝向矩阵，变换矩阵
//更改摄像机视角=更改物体位置和角度，操作为互逆矩阵
//摄像机变换是先旋转再平移，所以物体需要先平移后旋转，且都是逆矩阵
Matrix lookat(Vec3f eye, Vec3f center, Vec3f up) {
    //计算出z，根据z和up算出x，再算出y
    Vec3f z = (eye - center).normalize();
    Vec3f x = (up ^ z).normalize();
    Vec3f y = (z ^ x).normalize();
    Matrix rotation = Matrix::identity(4);
    Matrix translation = Matrix::identity(4);
    //***矩阵的第四列是用于平移的。因为观察位置从原点变为了center，所以需要将物体平移-center***
    for (int i = 0; i < 3; i++) {
        rotation[i][3] = -center[i];
    }
    //正交矩阵的逆 = 正交矩阵的转置
    //矩阵的第一行即是现在的x
    //矩阵的第二行即是现在的y
    //矩阵的第三行即是现在的z
    //***矩阵的三阶子矩阵是当前视线旋转矩阵的逆矩阵***
    for (int i = 0; i < 3; i++) {
        rotation[0][i] = x[i];
        rotation[1][i] = y[i];
        rotation[2][i] = z[i];
    }
    //这样乘法的效果是先平移物体，再旋转
    Matrix res = rotation*translation;
    return res;
}

//绘制三角形(坐标1，坐标2，坐标3，顶点光照强度1，顶点光照强度2，顶点光照强度3，tga指针，zbuffer)
void triangle(Vec3i t0, Vec3i t1, Vec3i t2, float ity0, float ity1, float ity2, Vec2i uv0, Vec2i uv1, Vec2i uv2,float dis0, float dis1, float dis2, TGAImage &image, int *zbuffer) {
    //按照y分割为两个三角形
    if (t0.y==t1.y && t0.y==t2.y) return;
    if (t0.y>t1.y) { std::swap(t0, t1); std::swap(ity0, ity1); std::swap(uv0, uv1);}
    if (t0.y > t2.y) { std::swap(t0, t2); std::swap(ity0, ity2); std::swap(uv0, uv2); }
    if (t1.y > t2.y) { std::swap(t1, t2); std::swap(ity1, ity2); std::swap(uv1, uv2); }
    int total_height = t2.y-t0.y;
    for (int i=0; i<total_height; i++) {
        bool second_half = i>t1.y-t0.y || t1.y==t0.y;
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height;
        //计算A,B两点的坐标
        Vec3i A    =               t0  + Vec3f(t2-t0  )*alpha;
        Vec3i B    = second_half ? t1  + Vec3f(t2-t1  )*beta : t0  + Vec3f(t1-t0  )*beta;
        //计算A,B两点的光照强度
        float ityA =               ity0 +   (ity2-ity0)*alpha;
        float ityB = second_half ? ity1 +   (ity2-ity1)*beta : ity0 +   (ity1-ity0)*beta;
        //计算UV
        Vec2i uvA = uv0 + (uv2 - uv0) * alpha;
        Vec2i uvB = second_half ? uv1 + (uv2 - uv1) * beta : uv0 + (uv1 - uv0) * beta;
        //计算距离
        float disA = dis0 + (dis2 - dis0) * alpha;
        float disB = second_half ? dis1 + (dis2 - dis1) * beta : dis0 + (dis1 - dis0) * beta;
        if (A.x>B.x) { std::swap(A, B); std::swap(ityA, ityB); }
        //x坐标作为循环控制
        for (int j=A.x; j<=B.x; j++) {
            float phi = B.x==A.x ? 1. : (float)(j-A.x)/(B.x-A.x);
            //计算当前需要绘制点P的坐标，光照强度
            Vec3i    P = Vec3f(A) +  Vec3f(B-A)*phi;
            float ityP =    ityA  + (ityB-ityA)*phi;
            ityP = std::min(1.f, std::abs(ityP)+0.01f);
            Vec2i uvP = uvA + (uvB - uvA) * phi;
            float disP = disA + (disB - disA) * phi;
            int idx = P.x+P.y*width;
            //边界限制
            if (P.x>=width||P.y>=height||P.x<0||P.y<0) continue;
            if (zbuffer[idx]<P.z) {
                zbuffer[idx] = P.z;
                TGAColor color = model->diffuse(uvP);
                image.set(P.x, P.y, TGAColor(color.bgra[2], color.bgra[1], color.bgra[0])*ityP*(20.f/std::pow(disP,2.f)));
                //image.set(P.x, P.y, TGAColor(255,255,255)* ityP);
            }
        }
    }
}
//计算某点的质心坐标,用来插值
Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
    Vec3f s[2];
    //计算[AB,AC,PA]的x和y分量，
    for (int i = 0; i < 2; i++) {
        s[i][0] = C[i] - A[i];//AC
        s[i][1] = B[i] - A[i];//AB
        s[i][2] = A[i] - P[i];//PA
    }
    //[u,v,1]和[AB,AC,PA]对应的x和y向量都垂直，所以叉乘,s[0]为x,s[1]为y
    Vec3f u = s[0] ^ s[1];
    //三点共线时，会导致u[2]为0，此时返回(-1,1,1)
    if (std::abs(u[2]) > 1e-2)//1e-2:0.01
        //若1-u-v，u，v全为大于0的数，表示点在三角形内部
        return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);//u[u,v,1]--->p[1-u-v,u,v];
    return Vec3f(-1, 1, 1);
}

bool insideTriangle(float x, float y, Vec3i t0, Vec3i t1, Vec3i t2) {
    Vec3f ab(t1.x - t0.x, t1.y - t0.y, 1);
    Vec3f bc(t2.x - t1.x, t2.y - t1.y, 1);
    Vec3f ca(t0.x - t2.x, t0.y - t2.y, 1);
    //Vec2f p0 = p-(Vec2f)t0;
    Vec3f ap(x - t0.x, y - t0.y, 1);
    Vec3f bp(x - t1.x, y - t1.y, 1);
    Vec3f cp(x - t2.x, y - t2.y, 1);

    if ((ab ^ (bc)) * (ab ^ (ap)) > 0 && (bc ^ (ca)) * (bc ^ (bp)) > 0 &&
        (ca ^ (ab)) * (ca ^ (cp)) > 0)
    {
        return true;
    }
    else
        return false;
}

void triangle_pixel(Vec3i t0, Vec3i t1, Vec3i t2, float ity0, float ity1, float ity2, Vec2i uv0, Vec2i uv1, Vec2i uv2, float dis0, float dis1, float dis2, TGAImage& image, int* zbuffer) {
    //找到bounding box,
    int min_x = std::min(t0.x, std::min(t1.x, t2.x));
    int min_y = std::min(t0.y, std::min(t1.y, t2.y));
    int max_x = std::max(t0.x, std::max(t1.x, t2.x));
    int max_y = std::max(t0.y, std::max(t1.y, t2.y));
    for (int x = min_x; x <= max_x; x++)
    {
        for (int y = min_y; y <= max_y; y++)
        {
            bool MSAA = false;
            Vec2f msaa[4] = {
                Vec2f(0.25, 0.25),
                Vec2f(0.75, 0.25),
                Vec2f(0.25, 0.75),
                Vec2f(0.75,0.75),
            };
            //MSAA实现不好，错误，所以先不用
            if (MSAA) {
                int cont = 0;
                for (int i = 0; i < 4; i++)
                {
                    if (insideTriangle(x + msaa[i].x, y + msaa[i].y, t0, t1, t2)) {
                        cont++;
                    }
                }
                if (cont > 0) {

                    //TGAColor colorNew(color.r * (cont / 4), color.g * (cont / 4), color.b * (cont / 4), color.a);
                   // image.set(x, y, colorNew);
                }
            }
            else {
                if (insideTriangle(x + 0.5, y + 0.5, t0, t1, t2)) {
                    Vec3f P;
                    P.x = x + 0.5;
                    P.y = y + 0.5;
                    Vec3f bc_coor = barycentric(t0, t1, t2, P);

                    float alpha = bc_coor.x;
                    float beta = bc_coor.y;
                    float gamma = bc_coor.z;

                    float inter_z = 1 / (alpha + gamma + beta);

                    float zx = alpha * t0.x + beta * t1.x + gamma * t2.x;
                    float zy = alpha * t0.y + beta * t1.y + gamma * t2.y;
                    float intensity_inter = (alpha * ity0 + beta * ity1 + gamma * ity2)*inter_z;
                    float distence_inter = (alpha * dis0 + beta * dis1 + gamma * dis2)*inter_z;
                    //注意此处，没有进行矫正，使用屏幕空间里的z值算的，若要矫正，需要世界坐标中的z值进行矫正！
                    float zp = alpha * t0.z + beta * t1.z + gamma * t2.z;

                    if (zbuffer[(int)(zx + zy * width)] < zp) {
                        zbuffer[(int)(zx + zy * width)] = zp;

                        Vec2i uvp = (uv0 * alpha + uv1 * beta + uv2 * gamma)*inter_z;
           
                        TGAColor color = model->diffuse(uvp);
                        //使用未校正过的alpha...值，使得矫正后代intensity和disp出现较大的误差
                        //image.set(x, y, TGAColor(color.r * intensity, color.g * intensity, color.b * intensity));
                        image.set(x, y, TGAColor(TGAColor(color.bgra[0], color.bgra[1], color.bgra[2]) *intensity_inter*  (20.f / std::pow(distence_inter, 2.f))));
                         //image.set(x, y, TGAColor(TGAColor(color.bgra[2], color.bgra[1], color.bgra[0]) ));
                    }
                }
            }
        }
    }

}

int main(int argc, char** argv) {
    //读取模型
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }
    //构造zbuffer并初始化
    zbuffer = new int[width*height];
    for (int i=0; i<width*height; i++) {
        zbuffer[i] = std::numeric_limits<int>::min();
    }
    //绘制模型
    {
        //模型变换矩阵
        Matrix ModelView  = lookat(eye, center, Vec3f(0,1,0));

        //透视矩阵
        Matrix Projection = Matrix::identity(4);
        Projection[3][2] = -1.f / (eye - center).norm();

        //视角矩阵
        Matrix ViewPort   = viewport(width/8, height/8, width*3/4, height*3/4);

        /*std::cerr << ModelView << std::endl;
        std::cerr << Projection << std::endl;
        std::cerr << ViewPort << std::endl;
        Matrix z = (ViewPort*Projection*ModelView);
        std::cerr << z << std::endl;*/

        TGAImage image(width, height, TGAImage::RGB);
        for (int i=0; i<model->nfaces(); i++) {
            std::vector<int> face = model->face(i);
            Vec3i screen_coords[3];
            float intensity[3];
            float distance[3];
            for (int j=0; j<3; j++) {
                Vec3f v = model->vert(face[j]);
                Matrix m_v = ModelView* Matrix(v);
                screen_coords[j] =  Vec3f(ViewPort*Projection* m_v);
                intensity[j] = model->norm(i, j)*light_dir;
                Vec3f new_v = Vec3f(m_v);
                distance[j] = std::pow((std::pow(new_v.x - eye.x,2.0f)+ std::pow(new_v.y - eye.y, 2.0f)+ std::pow(new_v.z - eye.z, 2.0f)),0.5f);
            }
            Vec2i uv[3];
            for (int k = 0; k < 3; k++) {
                uv[k] = model->uv(i, k);
            }
            //triangle(screen_coords[0], screen_coords[1], screen_coords[2], intensity[0], intensity[1], intensity[2], uv[0], uv[1], uv[2], distance[0], distance[1], distance[2], image, zbuffer);
            triangle_pixel(screen_coords[0], screen_coords[1], screen_coords[2], intensity[0], intensity[1], intensity[2], uv[0], uv[1], uv[2], distance[0], distance[1], distance[2], image, zbuffer);
        }
        image.flip_vertically();
        image.write_tga_file("new_output9.tga");
    }
    //输出zbuffer图像
    {
        TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                zbimage.set(i, j, TGAColor(zbuffer[i+j*width]));
            }
        }
        zbimage.flip_vertically();
        zbimage.write_tga_file("zbuffer.tga");
    }
    delete model;
    delete [] zbuffer;
    return 0;
}

