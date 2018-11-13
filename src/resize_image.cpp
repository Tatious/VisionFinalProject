#include <cmath>
#include "image.h"

using namespace std;

float Image::nn_interpolate(float x, float y, int c) const
  {
    return get_pixel(round(x), round(y), c);
  }

float Image::bilinear_interpolate(float x, float y, int c) const
  {
    int xL = ceil(x);
    int xS = floor(x);
    int yL = ceil(y);
    int yS = floor(y);

    float dx1 = x - xS;
    float dx2 = xL - x;
    float dy1 = y - yS;
    float dy2 = yL - y;

    float areaUL = dx1 * dy1;
    float areaUR = dx2 * dy1;
    float areaLL = dx1 * dy2;
    float areaLR = dx2 * dy2;

    float UL = get_pixel(xS, yS, c);
    float UR = get_pixel(xL, yS, c);
    float LL = get_pixel(xS, yL, c);
    float LR = get_pixel(xL, yL, c);

    if (xL == xS && yL == yS) {
      return UL;
    } else if (xL == xS) {
      return dy1 * LL + dy2 * UL;
    } else if (yL == yS) {
      return dx1 * UL + dx2 * UR;
    } else {
      return UL * areaLR + UR * areaLL + LL * areaUR + LR * areaUL;
    }

  }

float bilinear_interpolate(const Image& im, float x, float y, int c) { return im.bilinear_interpolate(x,y,c); }
float       nn_interpolate(const Image& im, float x, float y, int c) { return im.nn_interpolate(x,y,c); }

Image nn_resize(const Image& im, int w, int h)
  {
    if (w == 0 || h == 0)
      return make_image(1, 1, 1);

    Image ret = make_image(w,h,im.c);
    float xRatio = im.w / (float) w;
    float yRatio = im.h / (float) h;

    for (int k = 0; k < im.c; k++) {
      for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
          float x = (i + 0.5) * xRatio - 0.5;
          float y = (j + 0.5) * yRatio - 0.5;
          float resizePixel = nn_interpolate(im, x, y, k);
          set_pixel(ret, i, j, k, resizePixel);
        }
      }
    }
    return ret;
  }



Image bilinear_resize(const Image& im, int w, int h)
  {
    if (w == 0 || h == 0)
      return make_image(1, 1, 1);

    Image ret = make_image(w,h,im.c);
    float xRatio = im.w / (float) w;
    float yRatio = im.h / (float) h;

    for (int k = 0; k < im.c; k++) {
      for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
          float x = (i + 0.5) * xRatio - 0.5;
          float y = (j + 0.5) * yRatio - 0.5;
          float resizePixel = bilinear_interpolate(im, x, y, k);
          set_pixel(ret, i, j, k, resizePixel);
        }
      }
    }
    return ret;
  }
