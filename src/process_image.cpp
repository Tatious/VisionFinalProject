#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include "image.h"

using namespace std;

int clamp(int n, int min, int max);
Image rgb_to_grayscale(const Image& im)
  {
    assert(im.c == 3);
    Image gray = make_image(im.w, im.h, 1);
    float grayPixel;
    for (int j = 0; j < im.h; j++) {
      for (int i = 0; i < im.w; i++) {
        grayPixel = 0.299 * get_pixel(im, i, j, 0) + 0.587 * get_pixel(im, i, j, 1) + .114 * get_pixel(im, i, j, 2);
        set_pixel(gray, i, j, 0, grayPixel);
      }
    }
    return gray;
  }

void shift_image(Image& im, int c, float v)
  {
    for (int j = 0; j < im.h; j++) {
      for (int i = 0; i < im.w; i++) {
        float shiftPixel = get_pixel(im, i, j, c) + v;
        set_pixel(im, i, j, c, shiftPixel);
      }
    }
  }

void scale_image(Image& im, int c, float v)
  {
    for (int j = 0; j < im.h; j++) {
      for (int i = 0; i < im.w; i++) {
        float scalePixel = get_pixel(im, i, j, c) * v;
        set_pixel(im, i, j, c, scalePixel);
      }
    }
  }

void clamp_image(Image& im)
  {
    for (int k = 0; k < im.c; k++) {
      for (int j = 0; j < im.h; j++) {
        for (int i = 0; i < im.w; i++) {
          float clampPixel = get_pixel(im, i, j, k);
          clampPixel = clampPixel < 0 ? 0 : clampPixel > 1 ? 1 : clampPixel;
          set_pixel(im, i, j, k, clampPixel);
        }
      }
    }
  }

void Image::clamp(void) { clamp_image(*this); }
void Image::shift(int c, float v) { shift_image(*this,c,v); }
void Image::scale(int c, float v) { scale_image(*this,c,v); }


// These might be handy
float three_way_max(float a, float b, float c)
  {
  return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
  }

float three_way_min(float a, float b, float c)
  {
  return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
  }

void rgb_to_hsv(Image& im)
  {
    assert(im.c == 3);

    double R, G, B, V, C, S, H_prime, H;
    for (int j = 0; j < im.h; j++) {
      for (int i = 0; i < im.w; i++) {
        R = get_pixel(im, i, j, 0);
        G = get_pixel(im, i, j, 1);
        B = get_pixel(im, i, j, 2);

        V = three_way_max(R, G, B);
        C = (V - three_way_min(R, G, B));
        S = V == 0 ? 0 : C / V;
        H_prime = C == 0 ? 0 : V == R ? (G - B) / C : V == G ? ((B - R) / C) + 2 : ((R - G) / C) + 4;
        H = (H_prime / 6) + (H_prime < 0 ? 1 : 0);

        set_pixel(im, i, j, 0, H);
        set_pixel(im, i, j, 1, S);
        set_pixel(im, i, j, 2, V);
      }
    }
  }

void hsv_to_rgb(Image& im)
  {
    assert(im.c == 3);

    double H, S, V, C, X, R, G, B;
    for (int j = 0; j < im.h; j++) {
      for (int i = 0; i < im.w; i++) {
        H = get_pixel(im, i, j, 0);
        S = get_pixel(im, i, j, 1);
        V = get_pixel(im, i, j, 2);

        C = V * S;
        X = C * (1 - fabs(fmod(H * 6, 2) - 1));

        R = 0, G = 0, B = 0;
        switch ((int) (H * 6.0)) {
        case 0:
        R = C;
        G = X;
        break;
        case 1:
        R = X;
        G = C;
        break;
        case 2:
        G = C;
        B = X;
        break;
        case 3:
        G = X;
        B = C;
        break;
        case 4:
        R = X;
        B = C;
        break;
        default:
        R = C;
        B = X;
      }

      {
        double m = V - C;
        R += m;
        G += m;
        B += m;
      }

      set_pixel(im, i, j, 0, R);
      set_pixel(im, i, j, 1, G);
      set_pixel(im, i, j, 2, B);

      }
    }
  }


struct RGBcolor { float r,g,b; };
struct XYZcolor { float x,y,z; };
struct LCHcolor { float l,c,h; };

float l2g(float a)
  {
  if(a<0.0031308)return 12.92*a;
  else return 1.055*powf(a,1.0f/2.4f)-0.055;
  }

float g2l(float a)
  {
  if(a<0.04045)return a/12.92;
  else return powf((a+0.055f)/1.055f,2.4f);
  }

RGBcolor linear2gamma(RGBcolor a)
  {
  a.r=l2g(a.r);
  a.g=l2g(a.g);
  a.b=l2g(a.b);
  return a;
  }

RGBcolor gamma2linear(RGBcolor a)
  {
  a.r=g2l(a.r);
  a.g=g2l(a.g);
  a.b=g2l(a.b);
  return a;
  }

XYZcolor toXYZ(RGBcolor a)
  {
  XYZcolor b;
  a=gamma2linear(a);
  b.x=0.412383*a.r+0.357585*a.g+0.18048  *a.b;
  b.y=0.212635*a.r+0.71517 *a.g+0.072192 *a.b;
  b.z=0.01933 *a.r+0.119195*a.g+0.950528 *a.b;
  return b;
  }


RGBcolor toRGB(XYZcolor a)
  {
  RGBcolor b;
  b.r=(3.24103  )*a.x+(-1.53741 )*a.y +(-0.49862 )*a.z;
  b.g=(-0.969242)*a.x+(1.87596  )*a.y +(0.041555 )*a.z;
  b.b=(0.055632 )*a.x+(-0.203979)*a.y +(1.05698  )*a.z;
  b=linear2gamma(b);
  return b;
  }

LCHcolor rgb2lch(RGBcolor a)
  {
  LCHcolor b={0.f,0.f,0.f};
  XYZcolor c=toXYZ(a);

  if(c.x==0.f && c.y==0.f && c.z==0.f)return b;


  float u1=4*c.x/(1*c.x+15*c.y+3*c.z);
  float v1=9*c.y/(1*c.x+15*c.y+3*c.z);


  float un=0.2009;
  float vn=0.4610;

  float cutoff=powf(6.f/29.f,3);

  float l=0;
  if(c.y<=cutoff)l=powf(29.f/3.f,3)*c.y;
  else l=116.f*powf(c.y,1.f/3.f)-16.f;
  float u=13.f*l*(u1-un);
  float v=13.f*l*(v1-vn);


  b.l=l;
  b.c=sqrtf(u*u+v*v);
  b.h=atan2f(u,v);

  return b;
  }

RGBcolor lch2rgb(LCHcolor a)
  {
  XYZcolor b={0.f,0.f,0.f};

  if(a.l==0.f && a.c==0.f && a.h==0.f)return toRGB(b);

  float u=a.c*sinf(a.h);
  float v=a.c*cosf(a.h);
  float l=a.l;



  float un=0.2009;
  float vn=0.4610;

  float cutoff=8;


  float u1=u/(13.f*l)+un;
  float v1=v/(13.f*l)+vn;

  if(l<=cutoff)b.y=l*powf(3.f/29.f,3);
  else b.y=powf((l+16.f)/116.f,3);

  b.x=b.y*(9*u1)/(4*v1);
  b.z=b.y*(12-3*u1-20*v1)/(4*v1);

  //printf("xyz2   %f %f %f\n",b.x,b.y,b.z);

  return toRGB(b);
  }


void rgb_to_lch(Image& im)
  {
    assert(im.c == 3);

    float R, G, B;
    for (int j = 0; j < im.h; j++) {
      for (int i = 0; i < im.w; i++) {
        R = get_pixel(im, i, j, 0);
        G = get_pixel(im, i, j, 1);
        B = get_pixel(im, i, j, 2);
        RGBcolor rgb = {R, G, B};
        LCHcolor lch = rgb2lch(rgb);
        set_pixel(im, i, j, 0, lch.l);
        set_pixel(im, i, j, 1, lch.c);
        set_pixel(im, i, j, 2, lch.h);
      }
    }
  }

void lch_to_rgb(Image& im)
  {
    assert(im.c == 3);

    float L, C, H;
    for (int j = 0; j < im.h; j++) {
      for (int i = 0; i < im.w; i++) {
        L = get_pixel(im, i, j, 0);
        C = get_pixel(im, i, j, 1);
        H = get_pixel(im, i, j, 2);
        LCHcolor lch = {L, C, H};
        RGBcolor rgb = lch2rgb(lch);

        set_pixel(im, i, j, 0, rgb.r);
        set_pixel(im, i, j, 1, rgb.g);
        set_pixel(im, i, j, 2, rgb.b);
      }
    }
  }

void Image::HSVtoRGB(void) { hsv_to_rgb(*this); }
void Image::RGBtoHSV(void) { rgb_to_hsv(*this); }
void Image::LCHtoRGB(void) { lch_to_rgb(*this); }
void Image::RGBtoLCH(void) { rgb_to_lch(*this); }

int clamp(int n, int min, int max) {
  return n < min ? min : n >= max ? max - 1 : n;
}
