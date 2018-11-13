#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>

#include "image.h"
#include "matrix.h"


Image make_lucas_gx_filter() {
  Image ret = make_image(3, 3, 1);
  float gx[9] = {0, 0, 0, -1, 0, 1, 0, 0, 0};
  for (int i = 0; i < 9; i++) {
    ret.data[i] = gx[i];
  }
  return ret;
}


Image make_lucas_gy_filter() {
  Image ret = make_image(3, 3, 1);
  float gy[9] = {0, -1, 0, 0, 0, 0, 0, 1, 0};
  for (int i = 0; i < 9; i++) {
    ret.data[i] = gy[i];
  }
  return ret;
}


// Calculate the time-structure matrix of an Image pair.
// const Image& im: the input Image.
// const Image& prev: the previous Image in sequence.
// float s: sigma used for Gaussian smoothing the gradients
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          3rd channel is IxIy, 4th channel is IxIt, 5th channel is IyIt.
Image time_structure_matrix(const Image& im, const Image& prev, float s) {
  assert(im.c == 1 && prev.c == 1 && "Only for grayscale images");

  Image S(im.w, im.h, 5);

  Image ix_image = convolve_image(prev, make_lucas_gx_filter(), 0);
  Image iy_image = convolve_image(prev, make_lucas_gy_filter(), 0);
  Image it_image = sub_image(prev, im);

  for (int j = 0; j < im.h; j++) {
    for (int i = 0; i < im.w; i++) {
      float ix = ix_image(i, j);
      float iy = iy_image(i, j);
      float it = it_image(i, j);

      S(i, j, 0) = ix * ix;
      S(i, j, 1) = iy * iy;
      S(i, j, 2) = ix * iy;
      S(i, j, 3) = ix * it;
      S(i, j, 4) = iy * it;
    }
  }

  return fast_smooth_image(S, s);
}


// Compute the eigenvalues of the structure matrix
// Compute the eigenvalues only of S'S (the first three channels only)
// const Image& ts: the time-structure matrix
// returns: 2-channel image: 0-th channel : biggest eigenvalue,
//                           1-st channel : smallest
Image eigenvalue_matrix(const Image& ts) {

  Image im(ts.w,ts.h,2);

  for (int j = 0; j < im.h; j++) {
    for (int i = 0; i < im.w; i++) {
      float det = ts(i, j, 0) * ts(i, j, 1) - pow(ts(i, j, 2), 2);
      float tr = ts(i, j, 0) + ts(i, j, 1);

      float val = pow(tr, 2) / 4 - det;
      if (val < 0 && val > -1e-5) {
        val = 0.f;
      }

      float l1 = tr / 2.0 + pow(val, 0.5);
      float l2 = tr / 2.0 - pow(val, 0.5);

      im(i, j, 0) = max(l1, l2);
      im(i, j, 1) = min(l1, l2);
    }
  }

  return im;
}


vector<Image> make_image_pyramid(const Image& a, float factor, int levels) {
  assert(a.c==1 && "Only for grayscale");


  vector<Image> imgs(levels);
  imgs[0]=a;

  for(int l=1;l<levels;l++)
    {
    Image f=fast_smooth_image(imgs[l-1],factor/1.5);

    int bw=max(1,(int)(f.w/factor));
    int bh=max(1,(int)(f.h/factor));
    imgs[l]=bilinear_resize(f,bw,bh);
    }

  return imgs;
}


// Calculate the velocity given a structure Image
// const Image& S: time-structure Image
// const Image& ev: eigenvalue image
// Return: 2 channel (u,v) image  : the x and y
// velocities computed by inv(S'S)*(S'T)
Image velocity_image(const Image& S,const Image& ev) {
  Image v(S.w, S.h, 2);

  for (int j = 0; j < v.h; j++) {
    for (int i = 0; i < v.w; i++) {
      if (ev(i, j, 1) < 1e-5) {
        v(i, j, 0) = 0.0;
        v(i, j, 1) = 0.0;
      } else {
        Matrix2x2 StS = Matrix2x2(S(i, j, 0), S(i, j, 2), S(i, j, 2), S(i, j, 1));
        Matrix2x2 StS_inverted = StS.inverse();

        Vector2 StT = Vector2(S(i, j, 3), S(i, j, 4));
        Vector2 vect = StS_inverted * StT;

        v(i, j, 0) = vect.a;
        v(i, j, 1) = vect.b;
      }
    }
  }
  return v;
}


// Constrain the absolute value of each Image pixel
// const Image& im: Image to constrain
// float v: each pixel will be in range [-v, v]
void constrain_image(const Image& im, float v) {
  for(int i = 0; i < im.w*im.h*im.c; ++i)
    {
    if (im.data[i] < -v) im.data[i] = -v;
    if (im.data[i] >  v) im.data[i] =  v;
    }
}


// const Image& im: input image
// const Image& v: velocity image specifying how much each pixel moves
// return warped image, with same size as input image, as discussed on
// the github page.
Image warp_flow(const Image& im, const Image& v) {
  assert(im.c==1 && v.c==2 && "Only for grayscale and vel image needs 2 channels");
  assert(im.w==v.w && im.h==v.h && "Image and velocity need to be same size");

  float old_weight=1e-4;

  Image sum_weights(im.w, im.h, 1);
  Image sum_weighted_value(im.w, im.h, 1);

  for (int j = 0; j < im.h; j++) {
    for (int i = 0; i < im.w; i++) {
      sum_weights(i, j) = old_weight;
      sum_weighted_value(i, j) = old_weight * im(i, j);
    }
  }

  Image result(im.w, im.h, im.c);

  for (int j = 0; j < im.h; j++) {
    for (int i = 0; i < im.w; i++) {
      float x = i + v(i, j, 0);
      float y = j + v(i, j, 1);

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

      float orig = im(i, j);

      if (xL == xS && yL == yS) {
        set_pixel(sum_weights, x, y, 0, get_pixel(sum_weights, x, y, 0) + 1);
        set_pixel(sum_weighted_value, x, y, 0, get_pixel(sum_weighted_value, x, y, 0) + orig);
      } else if (xL == xS) {
        set_pixel(sum_weights, x, yS, 0, get_pixel(sum_weights, x, yS, 0) + dy2);
        set_pixel(sum_weights, x, yL, 0, get_pixel(sum_weights, x, yL, 0) + dy1);

        set_pixel(sum_weighted_value, x, yS, 0, get_pixel(sum_weighted_value, x, yS, 0) + orig * dy2);
        set_pixel(sum_weighted_value, x, yL, 0, get_pixel(sum_weighted_value, x, yL, 0) + orig * dy1);
      } else if (yL == yS) {
        set_pixel(sum_weights, xS, y, 0, get_pixel(sum_weights, xS, y, 0) + dx2);
        set_pixel(sum_weights, xL, y, 0, get_pixel(sum_weights, xL, y, 0) + dx1);

        set_pixel(sum_weighted_value, xS, y, 0, get_pixel(sum_weighted_value, xS, y, 0) + orig * dx2);
        set_pixel(sum_weighted_value, xL, y, 0, get_pixel(sum_weighted_value, xL, y, 0) + orig * dx1);
      } else {
        set_pixel(sum_weights, xL, yL, 0, get_pixel(sum_weights, xL, yL, 0) + areaUL);
        set_pixel(sum_weights, xL, yS, 0, get_pixel(sum_weights, xL, yS, 0) + areaLL);
        set_pixel(sum_weights, xS, yS, 0, get_pixel(sum_weights, xS, yS, 0) + areaLR);
        set_pixel(sum_weights, xS, yL, 0, get_pixel(sum_weights, xS, yL, 0) + areaUR);

        set_pixel(sum_weighted_value, xL, yL, 0, get_pixel(sum_weighted_value, xL, yL, 0) + orig * areaUL);
        set_pixel(sum_weighted_value, xL, yS, 0, get_pixel(sum_weighted_value, xL, yS, 0) + orig * areaLL);
        set_pixel(sum_weighted_value, xS, yS, 0, get_pixel(sum_weighted_value, xS, yS, 0) + orig * areaLR);
        set_pixel(sum_weighted_value, xS, yL, 0, get_pixel(sum_weighted_value, xS, yL, 0) + orig * areaUR);
      }
    }
  }

  for (int j = 0; j < im.h; j++) {
    for (int i = 0; i < im.w; i++) {
      result(i, j) = sum_weighted_value(i, j) / sum_weights(i, j);
    }
  }

  return result;
}


// Resize velocity image
// Image oldvel: old velocity
// int w,h : new sizes
// return new velocity image
Image velocity_resize(const Image& oldvel, int w, int h) {
  // TODO: resize the velocity image
  // do we just resize the image as normal or do we change the values as well?

  Image v = bilinear_resize(oldvel, w, h);
  v.scale(0, w / (float) oldvel.w);
  v.scale(1, h / (float) oldvel.h);
  return v;
};


void compute_iterative_pyramid_LK(LKIterPyramid& lk) {
  Image S;
  Image ev;
  Image v2;

  int h=lk.pyramid0[0].h;
  int w=lk.pyramid0[0].w;

  for(int q2=lk.pyramid_levels-1;q2>=0;q2--)
    {

    int pw=lk.pyramid1[q2].w;
    int ph=lk.pyramid1[q2].h;


    if(q2==lk.pyramid_levels-1)
      {
      lk.v=Image(pw,ph,2);
      lk.warped=lk.pyramid0[q2];
      }
    else
      {
      lk.v=velocity_resize(lk.v,pw,ph);
      lk.warped=warp_flow(lk.pyramid0[q2],lk.v);
      }

    for(int q1=0;q1<lk.lk_iterations;q1++)
      {
      S = time_structure_matrix(lk.pyramid1[q2], lk.warped, lk.smooth_structure);
      ev = eigenvalue_matrix(S);
      v2 = velocity_image(S, ev);

      v2=fast_smooth_image(v2,lk.smooth_vel);
      lk.v=lk.v+v2;

      constrain_image(lk.v,lk.clamp_vel);
      lk.warped=warp_flow(lk.pyramid0[q2],lk.v);
      }

    }



  lk.colorflow=vel2rgb(lk.v,lk.vel_color_scale);
  lk.error=(lk.warped-lk.pyramid1[0]).abs();

  if(lk.compute_all)
    {
    lk.all=Image(w*2,h*2,3);
    for(int c=0;c<3;c++)for(int q2=0;q2<h;q2++)for(int q1=0;q1<w;q1++)lk.all(q1+0,q2+0,c)=lk.t1(q1,q2,c);
    for(int c=0;c<3;c++)for(int q2=0;q2<h;q2++)for(int q1=0;q1<w;q1++)lk.all(q1+w,q2+0,c)=lk.colorflow(q1,q2,c);
    for(int c=0;c<3;c++)for(int q2=0;q2<h;q2++)for(int q1=0;q1<w;q1++)lk.all(q1+0,q2+h,c)=lk.warped(q1,q2);
    for(int c=0;c<3;c++)for(int q2=0;q2<h;q2++)for(int q1=0;q1<w;q1++)lk.all(q1+w,q2+h,c)=lk.error(q1,q2);
    }

  if(lk.compute_colored_ev)
    {
    lk.ev3=Image(ev.w,ev.h,3);
    memcpy(lk.ev3.data,ev.data,ev.size()*sizeof(float));
    }

}


// Calculate the optical flow between two images
// const Image& im: current Image
// Image prev: previous Image
// float smooth_win: amount to smooth structure matrix by
// float smooth_vel: amount to smooth velocity image
// returns: velocity matrix
Image optical_flow_images(const Image& im, const Image& prev, float smooth_win, float smooth_vel) {

  Image S = time_structure_matrix(im, prev, smooth_win);
  Image ev = eigenvalue_matrix(S);
  Image v = velocity_image(S, ev);
  if(smooth_vel==0)return v;
  return fast_smooth_image(v,smooth_vel);
}
