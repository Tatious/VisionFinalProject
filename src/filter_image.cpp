#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"

#define M_PI 3.14159265358979323846

void Image::l1_normalize(void) { ::l1_normalize(*this); }

void l1_normalize(Image& im)
  {
    float sum = 0;
    for (int k = 0; k < im.c; k++) {
      for (int j = 0; j < im.h; j++) {
        for (int i = 0; i < im.w; i++) {
          sum += get_pixel(im, i, j, k);
        }
      }
    }
    for (int k = 0; k < im.c; k++) {
      for (int j = 0; j < im.h; j++) {
        for (int i = 0; i < im.w; i++) {
          set_pixel(im, i, j, k, get_pixel(im, i, j, k) / sum);
        }
      }
    }
  }

Image make_box_filter(int w)
  {
    Image ret = make_image(w,w,1);
    for (int j = 0; j < ret.h; j++) {
      for (int i = 0; i < ret.w; i++) {
        set_pixel(ret, i, j, 0, 1);
      }
    }
    l1_normalize(ret);
    return ret;
  }

Image convolve_image(const Image& im, const Image& filter, int preserve)
  {
    assert(filter.c == im.c || filter.c == 1);

    Image ret = make_image(im.w, im.h, preserve == 1 ? im.c : 1);

    for (int j = 0; j < im.h; j++) {
      for (int i = 0; i < im.w; i++) {
        float sum = 0.0;
        for (int k = 0; k < im.c; k++) {
          for (int y = 0; y < filter.h; y++) {
            for (int x = 0; x < filter.w; x++) {
              int getX = i + (x - filter.w / 2);
              int getY = j + (y - filter.h / 2);

              float imageValue = get_pixel(im, getX, getY, k);
              float filterValue = get_pixel(filter, x, y, filter.c == im.c ? k : 0);

              sum += imageValue * filterValue;
            }
          }

          if (preserve == 1) {
            set_pixel(ret, i, j, k, sum);
            sum = 0.0;
          }
        }
        if (preserve != 1) {
          set_pixel(ret, i, j, 0, sum);
        }
      }
    }

    return ret;
  }

Image make_highpass_filter()
  {
    Image ret = make_image(3, 3, 1);
    float highpass[9] = {0, -1, 0, -1, 4, -1, 0, -1, 0};
    for (int i = 0; i < 9; i++) {
      ret.data[i] = highpass[i];
    }
    return ret;
  }

Image make_sharpen_filter()
  {
    Image ret = make_image(3, 3, 1);
    float sharpen[9] = {0, -1, 0, -1, 5, -1, 0, -1, 0};
    for (int i = 0; i < 9; i++) {
      ret.data[i] = sharpen[i];
    }
    return ret;
  }

Image make_emboss_filter()
  {
    Image ret = make_image(3, 3, 1);
    float emboss[9] = {-2, -1, 0, -1, 1, 1, 0, 1, 2};
    for (int i = 0; i < 9; i++) {
      ret.data[i] = emboss[i];
    }
    return ret;
  }

  // Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
  // Preserve sharpen & emboss since we want to keep the color channels in tact for the resulting image output.
  // For highpass, we want to flatten the image such that we get the resulting edges. We do not care about the color
  // in this case.

  // Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
  // Answer: Each of the above filters can produce values outside the range of 0 & 1, so we should clamp the values
  // before we try outputting any image.

Image make_gaussian_filter(float sigma)
  {
    int kernel_size = ceil(6 * sigma);
    kernel_size += kernel_size % 2 == 0 ? 1 : 0;

    Image ret = make_image(kernel_size, kernel_size, 1);

    float constant = 2 * M_PI * pow(sigma, 2);

    for (int x = 0; x < kernel_size; x++) {
      for (int y = 0; y < kernel_size; y++) {
        float kernX = (kernel_size / 2) - x;
        float kernY = (kernel_size / 2) - y;
        float val = exp(-(pow(kernX, 2) + pow(kernY, 2)) / (2 * pow(sigma, 2)));
        set_pixel(ret, x, y, 0, constant * val);
      }
    }
    l1_normalize(ret);
    return ret;
  }

Image add_image(const Image& a, const Image& b)
  {
    assert(a.w == b.w && a.h == b.h && a.c == b.c);

    Image ret = make_image(a.w, a.h, a.c);

    for (int k = 0; k < a.c; k++) {
      for (int j = 0; j < a.h; j++) {
        for (int i = 0; i < a.w; i++) {
          set_pixel(ret, i, j, k, get_pixel(a, i, j, k) + get_pixel(b, i, j, k));
        }
      }
    }
    return ret;
  }

Image sub_image(const Image& a, const Image& b)
  {
    assert(a.w == b.w && a.h == b.h && a.c == b.c);

    Image ret = make_image(a.w, a.h, a.c);

    for (int k = 0; k < a.c; k++) {
      for (int j = 0; j < a.h; j++) {
        for (int i = 0; i < a.w; i++) {
          set_pixel(ret, i, j, k, get_pixel(a, i, j, k) - get_pixel(b, i, j, k));
        }
      }
    }
    return ret;
  }

Image operator-(const Image& a, const Image& b) { return sub_image(a,b); }
Image operator+(const Image& a, const Image& b) { return add_image(a,b); }



Image make_gx_filter()
  {
    Image ret = make_image(3, 3, 1);
    float gx[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    for (int i = 0; i < 9; i++) {
      ret.data[i] = gx[i];
    }
    return ret;
  }

Image make_gy_filter()
  {
    Image ret = make_image(3, 3, 1);
    float gy[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    for (int i = 0; i < 9; i++) {
      ret.data[i] = gy[i];
    }
    return ret;
  }


void Image::feature_normalize(void) { ::feature_normalize(*this); }
void feature_normalize(Image& im)
  {
    for (int k = 0; k < im.c; k++) {
      float maxVal = -INFINITY;
      float minVal = INFINITY;

      for (int j = 0; j < im.h; j++) {
        for (int i = 0; i < im.w; i++) {
          maxVal = fmax(maxVal, get_pixel(im, i, j, k));
          minVal = fmin(minVal, get_pixel(im, i, j, k));
        }
      }

      float range = maxVal - minVal;

      for (int j = 0; j < im.h; j++) {
        for (int i = 0; i < im.w; i++) {
          float val = range == 0 ? 0 : (get_pixel(im, i, j, k) - minVal) / range;
          set_pixel(im, i, j, k, val);
        }
      }
    }


  }

void Image::feature_normalize_total(void) { ::feature_normalize_total(*this); }
void feature_normalize_total(Image& im)
  {
    float maxVal = -INFINITY;
    float minVal = INFINITY;

    for (int k = 0; k < im.c; k++) {
      for (int j = 0; j < im.h; j++) {
        for (int i = 0; i < im.w; i++) {
          maxVal = fmax(maxVal, get_pixel(im, i, j, k));
          minVal = fmin(minVal, get_pixel(im, i, j, k));
        }
      }
    }

    float range = maxVal - minVal;

    for (int k = 0; k < im.c; k++) {
      for (int j = 0; j < im.h; j++) {
        for (int i = 0; i < im.w; i++) {
          float val = range == 0 ? 0 : (get_pixel(im, i, j, k) - minVal) / range;
          set_pixel(im, i, j, k, val);
        }
      }
    }
  }

pair<Image,Image> sobel_image(const Image& im)
  {
    Image Mag(im.w,im.h);
    Image Theta(im.w,im.h);

    Image gxFilter = make_gx_filter();
    Image gyFilter = make_gy_filter();

    Image gx = convolve_image(im, gxFilter, 0);
    Image gy = convolve_image(im, gyFilter, 0);

    for (int y = 0; y < im.h; y++) {
      for (int x = 0; x < im.w; x++) {
        float gxPixel = get_pixel(gx, x, y, 0);
        float gyPixel = get_pixel(gy, x, y, 0);
        set_pixel(Mag, x, y, 0, sqrt(pow(gxPixel, 2) + pow(gyPixel, 2)));
        set_pixel(Theta, x, y, 0, atan2(gyPixel, gxPixel));
      }
    }

    return {Mag,Theta};
  }

Image colorize_sobel(const Image& im)
  {
    Image ret = make_image(im.w, im.h, im.c);
    std::pair<Image, Image> res = sobel_image(im);

    Image mag = res.first;
    Image theta = res.second;

    feature_normalize(mag);
    feature_normalize(theta);

    for (int j = 0; j < im.h; j++) {
      for (int i = 0; i < im.w; i++) {
        set_pixel(ret, i, j, 0, get_pixel(theta, i, j, 0));
        set_pixel(ret, i, j, 1, get_pixel(mag, i, j, 0));
        set_pixel(ret, i, j, 2, get_pixel(mag, i, j, 0));
      }
    }

    hsv_to_rgb(ret);

    return ret;
  }

Image bilateral_filter(const Image& im, float sigma1, float sigma2)
  {
    Image ret = make_image(im.w, im.h, im.c);
    float sigmaSpatial = 5.;
    float sigmaColor = 300.;
    float constant = 1. / (2. * M_PI * pow(sigmaSpatial, 2));

    int kernelSize = 7;
    Image filter = make_image(kernelSize, kernelSize, 1);

    for(int c = 0; c < im.c; c++){
      for(int y = 0; y < im.h; y++){
        for(int x = 0; x < im.w; x++){

          float sum = 0.;
          for (int i = -kernelSize / 2; i <= kernelSize / 2; i++) {
            for (int j = -kernelSize / 2; j <= kernelSize / 2; j++) {

              float spatial = constant * exp(-((pow(i, 2) + pow(j, 2)) / (2 * pow(sigmaSpatial, 2))));

              float colorDiff = get_pixel(im, x, y, c) - get_pixel(im, x + i, y + j, c);
              float color = (1. / (sigmaColor * sqrt(2 * M_PI))) * exp(-.5 * pow(colorDiff / sigmaColor, 2));

              sum += spatial * color;
              set_pixel(filter, i + kernelSize / 2, j + kernelSize / 2, 0, spatial * color);
            }
          }
          float finalValue = 0.;
          for (int i = -kernelSize / 2; i <= kernelSize / 2; i++) {
            for (int j = -kernelSize / 2; j <= kernelSize / 2; j++) {
              float originalValue = get_pixel(filter, i + kernelSize / 2, j + kernelSize / 2, 0);
              finalValue += get_pixel(im, x + i, y + j, c) * originalValue / sum;
            }
          }
          set_pixel(ret, x, y, c, finalValue);
        }
      }
    }

    return ret;
  }
