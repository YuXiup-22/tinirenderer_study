#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 800;

//�����㷨(����1������2��tgaָ�룬��ɫ)
void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
    bool steep = false;

    if (std::abs(p0.x-p1.x)<std::abs(p0.y-p1.y)) {
        std::swap(p0.x, p0.y);
        std::swap(p1.x, p1.y);
        steep = true;
    }
    if (p0.x>p1.x) {
        std::swap(p0, p1);
    }

    for (int x=p0.x; x<=p1.x; x++) {
        float t = (x-p0.x)/(float)(p1.x-p0.x);
        int y = p0.y*(1.-t) + p1.y*t;
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}

// Bresenham's line drawing algorithm
void Bresenham(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
    int dx, dy, x, y;
    //������ʼ��ֹ��,
    if (x0 > x1) {
        int temp = x1;
        x1 = x0;
        x0 = temp;
        temp = y1;
        y1 = y0;
        y0 = temp;
    }
    //ֱ��
    if (x0 == x1)
    {
        if (y0 < y1) {
            for (int i = y0; i <= y1; i++)
            {
                image.set(x0, i, color);
            }
        }
        else {
            for (int i = y1; i <= y0; i++)
            {
                image.set(x0, i, color);
            }
        }
    }
    if (y0 == y1) {
        for (int i = x0; i <= x1; i++)
        {
            image.set(i, y0, color);
        }
    }

    dx = x1 - x0;
    dy = y1 - y0;
    x = x0;
    y = y0;
    //trueΪб�ʵľ���ֵС��1
    bool m = (fabs(dx) > fabs(dy));
    //int dx1 = fabs(dx);
    //int dy1 = fabs(dy);
    int dx1 = dx;
    int dy1 = dy;
    int p1, p2, p;
    //����һ����
    image.set(x, y, color);
    if (m) {
        //б�ʴ���0С��1,��xΪ��
        if (y0 < y1) {
            p1 = 2 * dy1;
            p2 = 2 * dx1 - 2 * dy1;
            p = dx1 - 2 * dy1;
            for (int i = 0; i < dx1; i++)
            {
                if (p >= 0) {
                    p -= p1;
                    x++;
                }
                else {
                    p += p2;
                    x++;
                    y++;
                }
                image.set(x, y, color);
            }
        }
        //  -1<б��<0
        else {
            p = dx1 - 2 * dy1;
            p1 = 2 * (dx1 + dy1);
            p2 = 2 * dy1;
            for (int i = 0; i < dx1; i++)
            {
                if (p >= 0) {
                    p -= p1;
                    x++;
                    y--;
                }
                else {
                    p -= p2;
                    x++;
                }
                image.set(x, y, color);
            }
        }
    }
    else {
        //б�ʴ���1
        if (y0 < y1) {
            p = 2 * dx1 - dy1;
            p1 = 2 * (dx1 - dy1);
            p2 = 2 * dx1;
            for (int i = 0; i < dy1; i++)
            {
                if (p >= 0) {
                    p += p1;
                    x++;
                    y++;
                }
                else {
                    p += p2;
                    y++;
                }
                image.set(x, y, color);
            }
        }
        else {
            p = 2 * dx1 + dy1;
            p1 = 2 * (dx1 + dy1);
            p2 = 2 * dx1;
            for (int i = 0; i < -dy1; i++)
            {
                if (p >= 0) {
                    p += p1;
                    x++;
                    y--;
                }
                else {
                    p += p2;
                    y--;
                }
                image.set(x, y, color);
            }
        }
    }
}

//����������(����1������2������3��tgaָ�룬��ɫ)
void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
    //���������Ϊ0
    if (t0.y==t1.y && t0.y==t2.y) return;
    //����y�Ĵ�С�������������,�ǵ�t0��t2�ĸ߶ȵ���
    if (t0.y>t1.y) std::swap(t0, t1);
    if (t0.y>t2.y) std::swap(t0, t2);
    if (t1.y>t2.y) std::swap(t1, t2);
    int total_height = t2.y-t0.y;
    //�Ը߶Ȳ���Ϊѭ�����Ʊ�������ʱ����Ҫ����б�ʣ���Ϊ��ɫ���ÿ�ж��ᱻ���
    for (int i=0; i<total_height; i++) {
        //����t1�������ηָ�Ϊ����������
        bool second_half = i>t1.y-t0.y || t1.y==t0.y;//��true,����Ҫ����
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; 
        //����A,B���������
        Vec2i A =               t0 + (t2-t0)*alpha;
        Vec2i B = second_half ? t1 + (t2-t1)*beta : t0 + (t1-t0)*beta;
        if (A.x>B.x) std::swap(A, B);
        //����A,B�͵�ǰ�߶ȶ�tga��ɫ
        for (int j=A.x; j<=B.x; j++) {
            image.set(j, t0.y+i, color);
        }
    }
}


bool insideTriangle(float x, float y, Vec2i t0, Vec2i t1, Vec2i t2) {
    Vec3f ab(t1.x - t0.x, t1.y - t0.y, 1.0f);
    Vec3f bc(t2.x - t1.x, t2.y - t1.y, 1.0f);
    Vec3f ca(t0.x - t2.x, t0.y - t2.y, 1.0f);
    //Vec2f p0 = p-(Vec2f)t0;
    Vec3f ap(x - (float)t0.x, y - (float)t0.y,1.0f);
    Vec3f bp(x - (float)t1.x, y - (float)t1.y,1.0f);
    Vec3f cp(x - (float)t2.x, y - (float)t2.y,1.0f);

    if ((ab^(bc))*(ab^(ap)) > 0 && (bc^(ca))*(bc^(bp)) > 0 &&
        (ca^(ab))*(ca^(cp)) > 0)
    {
        return true;
    }else
        return false;
}
void triangle_pixel(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
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
                    if (insideTriangle((float)(x + msaa[i].x), (float)(y + msaa[i].y), t0, t1, t2)) {
                        cont++;
                    }
                }
                if (cont > 0) {
                    
                    TGAColor colorNew(color.r * (cont / 4), color.g * (cont / 4), color.b * (cont / 4), color.a);
                    image.set(x, y, colorNew);
                }
            }
            else {
            if (insideTriangle(x + 0.5, y + 0.5, t0, t1, t2)) {
                image.set(x, y, color);
            }
        }
        }
    }

}
int main(int argc, char** argv) {
    //��ȡmodel�ļ�
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }
    //����tgaimage
    TGAImage image(width, height, TGAImage::RGB);
   
    for (int i=0; i<model->nfaces(); i++) {
        //face��һ�����飬���ڴ洢һ�������������
        std::vector<int> face = model->face(i);
        Vec2i screen_coords[3];
        Vec3f world_coords[3];
        for (int j=0; j<3; j++) {//��ȡһ��face�е��������㣬�����䱾������������ NDC
            Vec3f v = model->vert(face[j]);
            //��Ļ����    (-1,-1)ӳ��Ϊ(0,0)  ��1,1��ӳ��Ϊ(width,height)
            screen_coords[j] = Vec2i((v.x+1.)*width/2., (v.y+1.)*height/2.);
            //��������    ��ģ������
            world_coords[j]  = v;
        }
        triangle_pixel(screen_coords[0], screen_coords[1], screen_coords[2], image, white);
        // //�������ɫ-FLAT SHADING����������������洹ֱ�������ǿ����1��ƽ����Ϊ0��
        // //ָ�����շ��򣬴�ʱdir�Ѿ������嵽��Դ�ķ���
   //     Vec3f light_dir(0,0,-1);
   //     //ģ�͵�����Ϊѭ�����Ʊ���
   ////     //������������㷨�����������������߲�˵���ķ���
   //     Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);
   //     //���߽��й�һ��
   //     n.normalize();
   ////     //ǿ��=������*���շ���   ���������͹��շ����غ�ʱ���������
   //     float intensity = n*light_dir;
   ////     //ǿ��С��0,�������ƽ����棬���䲻��
   //     if (intensity>0) {
   //         triangle_pixel(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity*255, intensity*255, intensity*255, 255));
   //     }
    }
    /*Vec2i t0[3] = { Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80) };
    Vec2i t1[3] = { Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180) };
    Vec2i t2[3] = { Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180) };
    triangle_pixel(t0[0], t0[1], t0[2], image, red);
    triangle_pixel(t1[0], t1[1], t1[2], image, white);
    triangle_pixel(t2[0], t2[1], t2[2], image, green);*/
    
    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output_MSAA2.tga");
    delete model;
    return 0;
}

