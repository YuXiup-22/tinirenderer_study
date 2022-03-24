#include "tgaimage.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    for (float t=0.; t<1.; t+=.1) {
        //t为步长，每次增加t*(x1-x0)或者t*(y1-y0)
        int x = x0*(1.-t) + x1*t;
        int y = y0*(1.-t) + y1*t;
        image.set(x, y, color);
    }
}
// Bresenham's line drawing algorithm
void Bresenham(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
    int dx, dy,x,y;
    //直线
    if (x0 == x1)
    {
        if (y0 < y1) 
        {
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
    //交换起始终止点
    if (x0 > x1) {
        int temp = x1;
        x1 = x0;
        x0 = temp;
        temp = y1;
        y1 = y0;
        y0 = temp;
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
            p1= 2 * (dx1 - dy1);
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
            p = -(2 * dx1 + dy1);
            p1 = -2 * (dx1 + dy1);
            p2 = -2 * dx1;
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
    TGAImage image(100, 100, TGAImage::RGB);
    //Bresenham(13, 20, 80, 40, image, white);
    //Bresenham(13, 80, 20, 20, image, white);
    //Bresenham(13, 20, 80, 40, image, white); //线段A
    //Bresenham(20, 13, 40, 80, image, red); //线段B
    //Bresenham(80, 40, 13, 20, image, red);//线段C
    //Bresenham(13, 80, 80, 60, image, red);
    //Bresenham(80, 10, 13, 80, image,red);
    //Bresenham(13, 80, 20, 10, image, red);
    Bresenham(20, 10, 13, 80, image, red);
    image.flip_vertically();
    image.write_tga_file("output6.tga");
    return 0;
}

