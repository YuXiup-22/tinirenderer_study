#include <vector>
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
Vec3f light_dir(0.2,0.15,-1);
Vec3f camera(0,0,3);
//�����λ��
Vec3f eye(2, 1, 3);
//����λ��
Vec3f center(0, 0, 1);
//������󣬱任����
//����������ӽ�=��������λ�úͽǶȣ�����Ϊ�������
//������任������ת��ƽ�ƣ�����������Ҫ��ƽ�ƺ���ת���Ҷ��������
Matrix lookat(Vec3f eye, Vec3f center, Vec3f up) {
    //�����z������z��up���x�������y
    Vec3f z = (eye - center).normalize();
    Vec3f x = (up ^ z).normalize();
    Vec3f y = (z ^ x).normalize();
    Matrix rotation = Matrix::identity(4);
    Matrix translation = Matrix::identity(4);
    //***����ĵ�����������ƽ�Ƶġ���Ϊ�۲�λ�ô�ԭ���Ϊ��center��������Ҫ������ƽ��-center***
    for (int i = 0; i < 3; i++) {
        rotation[i][3] = -center[i];
    }
    //����������� = ���������ת��
    //����ĵ�һ�м������ڵ�x
    //����ĵڶ��м������ڵ�y
    //����ĵ����м������ڵ�z
    //***����������Ӿ����ǵ�ǰ������ת����������***
    for (int i = 0; i < 3; i++) {
        rotation[0][i] = x[i];
        rotation[1][i] = y[i];
        rotation[2][i] = z[i];
    }
    //�����˷���Ч������ƽ�����壬����ת
    Matrix res = rotation * translation;
    return res;
}
//4d-->3d
//�������һ���������������һ������Ϊ0����ʾ������
//��Ϊ0����ʾ����
Vec3f m2v(Matrix m) {
    return Vec3f(m[0][0]/m[3][0], m[1][0]/m[3][0], m[2][0]/m[3][0]);
}

//3d-->4d
//���һ��1��ʾ����
Matrix v2m(Vec3f v) {
    Matrix m(4, 1);
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;
}

//�ӽǾ���
//������x��y����(-1,1)ת������Ļ����(100,700)    1/8width~7/8width
//zbuffer(-1,1)ת����0~255
Matrix viewport(int x, int y, int w, int h) {
    Matrix m = Matrix::identity(4);
    //��4�б�ʾƽ����Ϣ
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;
    //�Խ��߱�ʾ������Ϣ
    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
    return m;
}

//����������(�������꣬uv���꣬tgaָ�룬���ȣ�zbuffer)
void triangle(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, TGAImage &image, float intensity, int *zbuffer) {
    if (t0.y==t1.y && t0.y==t2.y) return;
    //�ָ������������
    if (t0.y>t1.y) { std::swap(t0, t1); std::swap(uv0, uv1); }
    if (t0.y>t2.y) { std::swap(t0, t2); std::swap(uv0, uv2); }
    if (t1.y>t2.y) { std::swap(t1, t2); std::swap(uv1, uv2); }
    //�ø߶���ѭ������
    int total_height = t2.y-t0.y;
    for (int i=0; i<total_height; i++) {
        //�ж�������һ������ȷ���߶�
        bool second_half = i>t1.y-t0.y || t1.y==t0.y;
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
        //���㵱ǰ�ı���
        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; // be careful: with above conditions no division by zero here
        //A��ʾt0��t2֮��ĵ�
        //B��ʾt0��t1֮��ĵ�
        Vec3i A   =               t0  + Vec3f(t2-t0  )*alpha;
        Vec3i B   = second_half ? t1  + Vec3f(t2-t1  )*beta : t0  + Vec3f(t1-t0  )*beta;
        //����UV
        
        Vec2i uvA =               uv0 +      (uv2-uv0)*alpha;
        Vec2i uvB = second_half ? uv1 +      (uv2-uv1)*beta : uv0 +      (uv1-uv0)*beta;
        //��֤B��A���ұ�
        if (A.x > B.x) { std::swap(A, B); }// std::swap(uvA, uvB);}
        //�ú�������Ϊѭ�����ƣ�����һ�н�����ɫ
        for (int j=A.x; j<=B.x; j++) {
            //���㵱ǰ����AB֮��ı���
            float phi = B.x==A.x ? 1. : (float)(j-A.x)/(float)(B.x-A.x);
            //�������ǰ�������,A��B������z����Ϣ
            Vec3i   P = Vec3f(A) + Vec3f(B-A)*phi;
            Vec2i uvP =     uvA +   (uvB-uvA)*phi;
            if (P.x < width && P.y < height)
            {
                //���㵱ǰzbuffer�±�=P.x+P.y*width
                int idx = P.x+P.y*width;
                //��ǰ���z����zbuffer��Ϣ�����ǵ���������zbuffer
                if (zbuffer[idx]<P.z) {
                    zbuffer[idx] = P.z;
                    TGAColor color = model->diffuse(uvP);
                    image.set(P.x, P.y, TGAColor(color.r*intensity, color.g*intensity, color.b*intensity));
                }
            }
        }
    }
}
//����ĳ�����������,������ֵ
Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
    Vec3f s[2];
    //����[AB,AC,PA]��x��y������
    for (int i = 0; i < 2; i++) {
        s[i][0] = C[i] - A[i];//AC
        s[i][1] = B[i] - A[i];//AB
        s[i][2] = A[i] - P[i];//PA
    }
    //[u,v,1]��[AB,AC,PA]��Ӧ��x��y��������ֱ�����Բ��,s[0]Ϊx,s[1]Ϊy
    Vec3f u = s[0]^s[1];
    //���㹲��ʱ���ᵼ��u[2]Ϊ0����ʱ����(-1,1,1)
    if (std::abs(u[2]) > 1e-2)//1e-2:0.01
        //��1-u-v��u��vȫΪ����0��������ʾ�����������ڲ�
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

void triangle_pixel(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, TGAImage& image,float intensity, int* zbuffer ) {
    //�ҵ�bounding box,
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
            //MSAAʵ�ֲ��ã����������Ȳ���
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
                    Vec3f bc_coor =barycentric(t0, t1, t2, P);
                    
                    float alpha = bc_coor.x;
                    float beta = bc_coor.y;
                    float gamma = bc_coor.z;
                    float zx = alpha * t0.x + beta * t1.x + gamma * t2.x;
                    float zy = alpha * t0.y + beta * t1.y + gamma * t2.y;

                    //ע��˴���û�н��н�����ʹ����Ļ�ռ����zֵ��ģ���Ҫ��������Ҫ���������е�zֵ���н�����
                    float zp = alpha * t0.z + beta * t1.z + gamma * t2.z;
                    
                    if (zbuffer[(int)(zx + zy * width)] < zp) {
                        zbuffer[(int)(zx + zy * width)] = zp;
                        Vec2i uvp = uv0 * alpha+uv1*beta+uv2*gamma;
                        TGAColor color = model->diffuse(uvp );
                        image.set(x, y, TGAColor(color.r * intensity, color.g * intensity, color.b * intensity));
                    }
                }
            }
        }
    }

}

int main(int argc, char** argv) {
    //��ȡģ��
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }
    //����zbuffer
    zbuffer = new int[width*height];
    for (int i=0; i<width*height; i++) {
        //��ʼ��zbuffer
        zbuffer[i] = std::numeric_limits<int>::min();
    }

    //����
    {
        ////ģ�ͱ任����
        //Matrix ModelView = lookat(eye, center, Vec3f(0, 1, 0));

        ////͸�Ӿ���
        //Matrix Projection = Matrix::identity(4);
        //Projection[3][2] = -1.f / (eye - center).norm();



        //��ʼ��ͶӰ����
        Matrix Projection = Matrix::identity(4);
        //��ʼ���ӽǾ���
        Matrix ViewPort   = viewport(width/8, height/8, width*3/4, height*3/4);
        ////ͶӰ����[3][2]=-1/c��cΪ���z����,�����е�����
        Projection[3][2] = -1.f/camera.z;
        //����tga
        TGAImage image(width, height, TGAImage::RGB);
        //��ģ����Ϊѭ�����Ʊ���
        for (int i=0; i<model->nfaces(); i++) {
            std::vector<int> face = model->face(i);
            Vec3i screen_coords[3];
            Vec3f world_coords[3];
            for (int j=0; j<3; j++) {
                //��ÿ����Ķ��㣬��ͶӰ
                Vec3f v = model->vert(face[j]);
                //�ӽǾ���*ͶӰ����*���꣬ͶӰ����ֱ�ӽ������ά������xoy�ϣ��ӽǾ�������ͶӰ����Ļ��
                screen_coords[j] = m2v(ViewPort * Projection * v2m(v));
                world_coords[j]  = v;
            }
            //���㷨����
            Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);
            n.normalize();
            //�������
            float intensity = n*light_dir;
            intensity = std::min(std::abs(intensity),1.f);
            if (intensity>0) {
                Vec2i uv[3];
                for (int k=0; k<3; k++) {
                    //uv(i,k)----��i����ĵ�k�������uvֵ
                    uv[k] = model->uv(i, k);
                }
                //����������
                triangle_pixel(screen_coords[0], screen_coords[1], screen_coords[2], uv[0], uv[1], uv[2], image, intensity, zbuffer);
               // triangle(screen_coords[0], screen_coords[1], screen_coords[2], uv[0], uv[1], uv[2], image, intensity, zbuffer);
            }
        }
        //tgaĬ��ԭ�������ϣ��ָ�Ϊ����
        image.flip_vertically();
        image.write_tga_file("new_output5.tga");
    }
    //���zbuffer
    {
        TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
        for (int i=0; i<width; i++) {
            for (int j=0; j<height; j++) {
                //���Խǳ��Խ��
                zbimage.set(i, j, TGAColor(zbuffer[i+j*width], 1));
            }
        }
        zbimage.flip_vertically();
        zbimage.write_tga_file("zbuffer.tga");
    }
    delete model;
    delete [] zbuffer;
    return 0;
}

