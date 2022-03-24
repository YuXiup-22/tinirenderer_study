#include <vector>
#include <cmath>
#include "tgaimage.h"   //tga画图库
#include "model.h"      //模型类，主要实现模型的读取
#include "geometry.h"   //几何库，主要定义了Vec2和Vec3类型

//定义颜色
const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
Model *model = NULL;
//定义宽度高度
const int width  = 800;
const int height = 800;

//画线算法(坐标1，坐标2，tga目标指针，指定颜色)
void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    //判断线段斜率的绝对值是否大于1
    bool steep = false;
    //大于1置为true，交换坐标各自的x和y。即变换为关于y=x或y=-x对称的点
    if (std::abs(x0-x1)<std::abs(y0-y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    //保证坐标2的x,y大于坐标1的x,y。
    if (x0>x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    //此时x1大于x0，且斜率在-1到1之间，用x做循环控制变量
    for (int x=x0; x<=x1; x++) {
        //根据x计算线段对应的y
        float t = (x-x0)/(float)(x1-x0);
        int y = y0*(1.-t) + y1*t;
        //若斜率大于1，真实坐标为(y,x)；否则为(x,y)
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
    //交换起始终止点,
    if (x0 > x1) {
        int temp = x1;
        x1 = x0;
        x0 = temp;
        temp = y1;
        y1 = y0;
        y0 = temp;
    }
    //直线
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
        for (int i = x0; i <=x1; i++)
        {
            image.set(i, y0, color);
        }
    }
    
    dx = x1 - x0;
    dy = y1 - y0;
    x = x0;
    y = y0;
    //true为斜率的绝对值小于1
    bool m = (fabs(dx) > fabs(dy));
    //int dx1 = fabs(dx);
    //int dy1 = fabs(dy);
    int dx1 = dx;
    int dy1 = dy;
    int p1, p2, p;
    //画第一个点
    image.set(x, y, color);
    if (m) {
        //斜率大于0小于1,以x为主
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
        //  -1<斜率<0
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
        //斜率大于1
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


int main(int argc, char** argv) {
    //命令行控制方式和代码方式构造model
    //构造模型(obj文件路径)
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
        //model = new Model("obj/cube.obj");
    }
    //构造tga(宽，高，指定颜色空间)  屏幕空间
    TGAImage image(width, height, TGAImage::RGB);
    //模型的面作为循环控制变量
    
    for (int i=0; i<model->nfaces(); i++) {
        //创建face数组用于保存一个face的三个顶点坐标
        std::vector<int> face = model->face(i);
        //遍历一个face的三个顶点，并两两连接，最后构成三角形
        for (int j=0; j<3; j++) {
            //顶点v0
            Vec3f v0 = model->vert(face[j]);
            //顶点v1
            Vec3f v1 = model->vert(face[(j+1)%3]);
            //根据顶点v0和v1画线
            //先要进行模型坐标到屏幕坐标的转换，就是NDC到屏幕坐标转换
            //(-1,-1)对应(0,0)   (1,1)对应(width,height)
            int x0 = (v0.x+1.)*width/2.;
            int y0 = (v0.y+1.)*height/2.;
            int x1 = (v1.x+1.)*width/2.;
            int y1 = (v1.y+1.)*height/2.;
            //画线
           // line(x0,y0, x1,y1, image, white);
            Bresenham(x0, y0, x1, y1, image, white);
        }
    }

    //tga默认原点在左上角，现需要指定为左下角，所以进行竖直翻转
    image.flip_vertically();
    image.write_tga_file("Bre_African_output1.tga");
    delete model;
    return 0;
}

