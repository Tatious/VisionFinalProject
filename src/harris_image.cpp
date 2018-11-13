#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>

#include "image.h"
//#include "matrix.h"

using namespace std;


// Create a feature descriptor for an index in an image.
// const Image& im: source image.
// int x,y: coordinates for the pixel we want to describe.
// returns: Descriptor for that index.
Descriptor describe_index(const Image& im, int x, int y, int w)
  {
  Descriptor d;
  d.p={(double)x,(double)y};
  d.data.reserve(w*w*im.c);

  // If you want you can experiment with other descriptors
  // This subtracts the central value from neighbors
  // to compensate some for exposure/lighting changes.
  for(int c=0;c<im.c;c++)
    {
    float cval = im.get_pixel(x,y,c);
    for(int dx=-w/2;dx<=w/2;dx++)for(int dy=-w/2;dy<=w/2;dy++)
      d.data.push_back(im.get_pixel(x+dx,y+dy,c)-cval);
    }
  return d;
  }

// Marks the spot of a point in an image.
// Image& im: image to mark.
// Point p: spot to mark in the image.
void mark_spot(Image& im, const Point& p)
  {
  int x = p.x;
  int y = p.y;

  for(int i = -9; i < 10; ++i)
    {
    im.set_pixel(x+i, y, 0, 1);
    im.set_pixel(x, y+i, 0, 1);
    im.set_pixel(x+i, y, 1, 0);
    im.set_pixel(x, y+i, 1, 0);
    im.set_pixel(x+i, y, 2, 1);
    im.set_pixel(x, y+i, 2, 1);
    }
  }

// Marks corners denoted by an array of descriptors.
// Image& im: image to mark.
// vector<Descriptor> d: corners in the image.
Image mark_corners(const Image& im, const vector<Descriptor>& d)
  {
  Image im2=im;
  for(auto&e1:d)mark_spot(im2,e1.p);
  return im2;
  }

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row Image of the filter.
Image make_1d_gaussian(float sigma)
  {
    //TIME(1);
    int kernel_size = ceil(6 * sigma);
    kernel_size += kernel_size % 2 == 0 ? 1 : 0;

    Image ret = make_image(kernel_size, 1, 1);
    double constant = (1. / (sigma * sqrt(2 * M_PI)));

    for (int i = -kernel_size / 2; i <= kernel_size / 2; i++) {
      set_pixel(ret, i + kernel_size / 2, 0, 0, constant * exp(- .5 * pow(i / sigma, 2)));
    }
    l1_normalize(ret);

    return ret;
  }

// Smooths an image using separable Gaussian filter.
// const Image& im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed Image.
Image smooth_image(const Image& im, float sigma)
  {
    //TIME(1);
    Image gx = make_1d_gaussian(sigma);
    Image gy(1, gx.w);
    for (int i = 0; i < gx.w; i++) {
      gy(0, i) = gx(i, 0);
    }
    Image ret = convolve_image(im, gx, 1);
    ret = convolve_image(ret, gy, 1);

    return ret;
  }

// Calculate the structure matrix of an image.
// const Image& im im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
Image structure_matrix(const Image& im2, float sigma)
  {
  //TIME(1);
  // only grayscale or rgb
  assert((im2.c==1 || im2.c==3) && "only grayscale or rgb supported");
  Image im;
  // convert to grayscale if necessary
  if(im2.c==1)im=im2;
  else im=rgb_to_grayscale(im2);

  Image S(im.w, im.h, 3);

  Image ix_image = convolve_image(im, make_gx_filter(), 0);
  Image iy_image = convolve_image(im, make_gy_filter(), 0);

  for (int j = 0; j < im.h; j++) {
    for (int i = 0; i < im.w; i++) {
      float ix = ix_image(i, j);
      float iy = iy_image(i, j);

      S(i, j, 0) = ix * ix;
      S(i, j, 1) = iy * iy;
      S(i, j, 2) = ix * iy;
    }
  }

  S = smooth_image(S, sigma);
  return S;
  }

// Estimate the cornerness of each pixel given a structure matrix S.
// const Image& im S: structure matrix for an image.
// returns: a response map of cornerness calculations.
// int method: 0: det(S)/tr(S)    1 (optional) : exact 2nd eigenvalue
Image cornerness_response(const Image& S, int method)
  {
  Image R(S.w, S.h);

  for (int j = 0; j < S.h; j++) {
    for (int i = 0; i < S.w; i++) {
      float det = S(i, j, 0) * S(i, j, 1) - pow(S(i, j, 2), 2);
      float tr = S(i, j, 0) + S(i, j, 1);
      if (method == 0) {
        R(i, j) = det / tr;
      } else if (method == 1) {
        R(i, j) = tr / 2.0 - pow((pow(tr, 2) / 4 - det), 0.5);
      }
    }
  }
  return R;
  }

// Perform non-max supression on an image of feature responses.
// const Image& im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: Image with only local-maxima responses within w pixels.
Image nms_image(const Image& im, int w)
  {
  //TIME(1);
  Image r = im;

  w = 2 * w + 1;
  for (int j = 0; j < im.h; j++) {
    for (int i = 0; i < im.w; i++) {

      bool non_max = false;
      for (int wx = -w / 2; wx <= w / 2 && !non_max; wx++) {
        for (int wy = -w / 2; wy <= w / 2 && !non_max; wy++) {
          if (im(i,j) < get_pixel(im, i + wx, j + wy, 0)) {
            non_max = true;
          }
        }
      }
      r(i, j) = non_max ? -999999 : im(i, j);
    }
  }
  return r;
  }


// Perform corner detection and extract features from the corners.
// const Image& im: input image.
// const Image& nms: nms image
// float thresh: threshold for cornerness.
// returns: vector of descriptors of the corners in the image.
vector<Descriptor> detect_corners(const Image& im, const Image& nms, float thresh, int window)
  {
  vector<Descriptor> d;

  for (int j = 0; j < im.h; j++) {
    for (int i = 0; i < im.w; i++) {
      if (nms(i, j) >= thresh) {
        d.push_back(describe_index(im, i, j, window));
      }
    }
  }
  return d;
  }

// Perform harris corner detection and extract features from the corners.
// const Image& im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// returns: vector of descriptors of the corners in the image.
vector<Descriptor> harris_corner_detector(const Image& im, float sigma, float thresh, int window, int nms, int corner_method)
  {
  // Calculate structure matrix
  Image S = structure_matrix(im, sigma);

  // Estimate cornerness
  Image R = cornerness_response(S,corner_method);

  // Run NMS on the responses
  Image Rnms = nms_image(R, nms);

  return detect_corners(im, Rnms, thresh, window);
  }

// Find and draw corners on an image.
// Image& im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
Image detect_and_draw_corners(const Image& im, float sigma, float thresh, int window, int nms, int corner_method)
  {
  vector<Descriptor> d = harris_corner_detector(im, sigma, thresh, window, nms, corner_method);
  return mark_corners(im, d);
  }
