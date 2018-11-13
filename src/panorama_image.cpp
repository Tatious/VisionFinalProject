#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>

#include "image.h"
#include "matrix.h"

#include <set>

using namespace std;

// Place two images side by side on canvas, for drawing matching pixels.
// const Image& a, b: images to place.
// returns: image with both a and b side-by-side.
Image both_images(const Image& a, const Image& b)
  {
  assert(a.c==b.c);
  Image both(a.w + b.w, a.h > b.h ? a.h : b.h, a.c);

  for(int k = 0; k < both.c; ++k)
    for(int j = 0; j < a.h; ++j)
      for(int i = 0; i < a.w; ++i)
        both(i, j, k) = a(i, j, k);

  for(int k = 0; k < both.c; ++k)
    for(int j = 0; j < b.h; ++j)
      for(int i = 0; i < b.w; ++i)
        both(i+a.w, j, k) = b(i, j, k);
  return both;
  }

// Draws lines between matching pixels in two images.
// const Image& a, b: two images that have matches.
// const vector<Match>& matches: array of matches between a and b.
// int inliers: number of inliers at beginning of matches, drawn in green.
// returns: image with matches drawn between a and b on same canvas.
Image draw_matches(const Image& a, const Image& b, const vector<Match>& matches, const vector<Match>& inliers)
  {
  Image both = both_images(a, b);

  for(int i = 0; i < (int)matches.size(); ++i)
    {
    int bx = matches[i].a->p.x;
    int ex = matches[i].b->p.x;
    int by = matches[i].a->p.y;
    int ey = matches[i].b->p.y;
    for(int j = bx; j < ex + a.w; ++j)
      {
      int r = (float)(j-bx)/(ex+a.w - bx)*(ey - by) + by;
      both.set_pixel(j, r, 0, 1);
      both.set_pixel(j, r, 1, 0);
      both.set_pixel(j, r, 2, 0);
      }
    }
  for(int i = 0; i < (int)inliers.size(); ++i)
    {
    int bx = inliers[i].a->p.x;
    int ex = inliers[i].b->p.x;
    int by = inliers[i].a->p.y;
    int ey = inliers[i].b->p.y;
    for(int j = bx; j < ex + a.w; ++j)
      {
      int r = (float)(j-bx)/(ex+a.w - bx)*(ey - by) + by;
      both.set_pixel(j, r, 0, 0);
      both.set_pixel(j, r, 1, 1);
      both.set_pixel(j, r, 2, 0);
      }
    }
  return both;
  }

// Draw the matches with inliers in green between two images.
// const Image& a, b: two images to match.
// vector<Match> m: matches
// Matrix H: the current homography
// thresh: for thresholding inliers
Image draw_inliers(const Image& a, const Image& b, const Matrix& H, const vector<Match>& m, float thresh)
  {
  vector<Match> inliers = model_inliers(H, m, thresh);
  Image lines = draw_matches(a, b, m, inliers);
  return lines;
  }

// Find corners, match them, and draw them between two images.
// const Image& a, b: images to match.
// float sigma: gaussian for harris corner detector. Typical: 2
// float thresh: threshold for corner/no corner. Typical: 1-5
// int nms: window to perform nms on. Typical: 3
Image find_and_draw_matches(const Image& a, const Image& b, float sigma, float thresh, int window, int nms, int corner_method)
  {
  vector<Descriptor> ad= harris_corner_detector(a, sigma, thresh, window, nms, corner_method);
  vector<Descriptor> bd= harris_corner_detector(b, sigma, thresh, window, nms, corner_method);
  vector<Match> m = match_descriptors(ad, bd);


  Image A=mark_corners(a, ad);
  Image B=mark_corners(b, bd);
  Image lines = draw_matches(A, B, m, {});

  return lines;
  }

// Calculates L1 distance between two floating point arrays.
// vector<float>& a,b: arrays to compare.
// returns: l1 distance between arrays (sum of absolute differences).
float l1_distance(const vector<float>& a,const vector<float>& b)
  {
  assert(a.size() == b.size() && "Arrays must have same size\n");

  float distance = 0.0;

  for (int i = 0; i < a.size(); i++) {
    distance += fabs(a[i] - b[i]);
  }

  return distance;
  }


// Finds best matches between descriptors of two images.
// const vector<Descriptor>& a, b: array of descriptors for pixels in two images.
// returns: best matches found. For each element in a[] find the index of best match in b[]
vector<int> match_descriptors_a2b(const vector<Descriptor>& a, const vector<Descriptor>& b)
  {
  vector<int> ind;

  for (int j = 0; j < (int) a.size(); j++) {
    int bind = -1; // <- find the best match (-1: no match)
    float best_distance = 1e10f;  // <- best distance

    for (int i = 0; i < (int) b.size(); i++) {
      float distance = l1_distance(a[j].data, b[i].data);
      if (distance < best_distance) {
        best_distance = distance;
        bind = i;
      }
    }

    ind.push_back(bind);
  }

  return ind;
  }



// Finds best matches between descriptors of two images.
// const vector<Descriptor>& a, b: array of descriptors for pixels in two images.
// returns: best matches found. each descriptor in a should match with at most
//          one other descriptor in b.
vector<Match> match_descriptors(const vector<Descriptor>& a, const vector<Descriptor>& b)
  {
  if (a.size() == 0 || b.size() == 0) return {};

  vector<Match> m;

  vector<int> a2b = match_descriptors_a2b(a,b);
  vector<int> b2a = match_descriptors_a2b(b,a);

  for (int i = 0; i < a2b.size(); i++) {
    if (i == b2a[a2b[i]]) {
      Match match = {&a[i], &b[a2b[i]], l1_distance(a[i].data, b[a2b[i]].data)};
      m.push_back(match);
    }
  }

  return m;
  }


// Apply a projective transformation to a point.
// const Matrix& H: homography to project point.
// const Point& p: point to project.
// returns: point projected using the homography.
Point project_point(const Matrix& H, const Point& p)
  {
  assert(H.rows == 3 && H.cols == 3);

  Matrix c(3,1);
  double* data = new double [3]{p.x, p.y, 1.0};
  c.data = data;

  Matrix result = H * c;

  return Point(result.data[0] / result.data[2],
                  result.data[1] / result.data[2]);
  }

// Calculate L2 distance between two points.
// const Point& p, q: points.
// returns: L2 distance between them.
double point_distance(const Point& p, const Point& q)
  {
  return pow(pow(p.x - q.x, 2) + pow(p.y - q.y, 2), 0.5);
  }

// Count number of inliers in a set of matches. Should also bring inliers
// to the front of the array.
// const Matrix& H: homography between coordinate systems.
// const vector<Match>& m: matches to compute inlier/outlier.
// float thresh: threshold to be an inlier.
// returns: inliers whose projected point falls within thresh of their match in the other image.
vector<Match> model_inliers(const Matrix& H, const vector<Match>& m, float thresh)
  {
  vector<Match> inliers;
  for (int i = 0; i < m.size(); i++) {
    Point transformed = project_point(H, m[i].a->p);
    if (point_distance(transformed, m[i].b->p) < thresh) {
      inliers.push_back(m[i]);
    }
  }

  return inliers;
  }

// Randomly shuffle matches for RANSAC.
// vector<Match>& m: matches to shuffle in place.
void randomize_matches(vector<Match>& m)
  {
    for (int i = m.size() - 1; i > 0; i--) {
      int j = rand() % (i + 1);
      Match temp = m[j];
      m[j] = m[i];
      m[i] = temp;
    }
  }

// Computes homography between two images given matching pixels.
// const vector<Match>& matches: matching points between images.
// int n: number of matches to use in calculating homography.
// returns: matrix representing homography H that maps image a to image b.
Matrix compute_homography_ba(const vector<Match>& matches) {
  if (matches.size() < 4) printf("Need at least 4 points for homography! %zu supplied\n",matches.size());
  if (matches.size() < 4) return Matrix::identity(3, 3);

  Matrix M(matches.size() * 2, 8);
  Matrix b(matches.size() * 2);

  M.data = new double [matches.size() * 2 * 8];
  b.data = new double [matches.size() * 2];

  for(int i = 0; i < (int)matches.size(); ++i) {
    double mx = matches[i].a->p.x;
    double my = matches[i].a->p.y;

    double nx = matches[i].b->p.x;
    double ny = matches[i].b->p.y;

    double val1[8] = {mx, my, 1, 0, 0, 0, -mx * nx, -my * nx};
    std::copy(val1, val1 + 8, M.data + (8 * (2 * i)));
    double val2[8] = {0, 0, 0, mx, my, 1, -mx * ny, -my * ny};
    std::copy(val2, val2 + 8, M.data + (8 * (2 * i + 1)));

    b.data[i * 2] = nx;
    b.data[i * 2 + 1] = ny;
  }

  Matrix a = solve_system(M, b);

  Matrix Hba(3, 3);
  Hba.data = new double [9];
  std::copy(a.data, a.data + 8, Hba.data);
  Hba.data[8] = 1.0;
  return Hba;
}

// Perform RANdom SAmple Consensus to calculate homography for noisy matches.
// vector<Match> m: set of matches.
// float thresh: inlier/outlier distance threshold.
// int k: number of iterations to run.
// int cutoff: inlier cutoff to exit early.
// returns: matrix representing most common homography between matches.
Matrix RANSAC(vector<Match> m, float thresh, int k, int cutoff) {
  if (m.size() < 4) printf("Need at least 4 points for RANSAC! %zu supplied\n",m.size());
  if (m.size() < 4) return Matrix::identity(3,3);
  //TIME(1);

  int best = 0;
  Matrix Hba = Matrix::translation_homography(256, 0);

  for (int i = 0; i < k; i++) {
    randomize_matches(m);
    vector<Match> homography_matches;
    for(int i = 0; i < 4; i++) { // (how many??)
      homography_matches.push_back(m[i]);
    }

    Matrix homography = compute_homography_ba(homography_matches);

    vector<Match> inliers = model_inliers(homography, m, thresh);
    if (inliers.size() > best) {
      best = inliers.size();
      Hba = compute_homography_ba(inliers);
      if (best > cutoff) {
        return Hba;
      }
    }

  }
  return Hba;
}

Image trim_image(const Image& a)
  {
  int minx=a.w-1;
  int maxx=0;
  int miny=a.h-1;
  int maxy=0;

  for(int q3=0;q3<a.c;q3++)for(int q2=0;q2<a.h;q2++)for(int q1=0;q1<a.w;q1++)if(a(q1,q2,q3))
    {
    minx=min(minx,q1);
    maxx=max(maxx,q1);
    miny=min(miny,q2);
    maxy=max(maxy,q2);
    }

  if(maxx<minx || maxy<miny)return a;

  Image b(maxx-minx+1,maxy-miny+1,a.c);

  for(int q3=0;q3<a.c;q3++)for(int q2=miny;q2<=maxy;q2++)for(int q1=minx;q1<=maxx;q1++)
    b(q1-minx,q2-miny,q3)=a(q1,q2,q3);

  return b;
  }

// MIN MAX MACROS
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


// Stitches two images together using a projective transformation.
// const Image& a, b: images to stitch.
// Matrix H: homography from image a coordinates to image b coordinates.
// float acoeff: blending coefficient
// returns: combined image stitched together.
Image combine_images(const Image& a, const Image& b, const Matrix& Hba, float ablendcoeff) {
  //TIME(1);

  Matrix Hinv=Hba.inverse();

  // Project the corners of image b into image a coordinates.
  Point c1 = project_point(Hinv, Point(0,0));
  Point c2 = project_point(Hinv, Point(b.w-1, 0));
  Point c3 = project_point(Hinv, Point(0, b.h-1));
  Point c4 = project_point(Hinv, Point(b.w-1, b.h-1));

  // Find top left and bottom right corners of image b warped into image a.
  Point topleft, botright;
  botright.x = max(c1.x, max(c2.x, max(c3.x, c4.x)));
  botright.y = max(c1.y, max(c2.y, max(c3.y, c4.y)));
  topleft.x = min(c1.x, min(c2.x, min(c3.x, c4.x)));
  topleft.y = min(c1.y, min(c2.y, min(c3.y, c4.y)));

  // Find how big our new image should be and the offsets from image a.
  int dx = min(0, (int)topleft.x);
  int dy = min(0, (int)topleft.y);
  int w = max(a.w, (int)botright.x) - dx;
  int h = max(a.h, (int)botright.y) - dy;

  //printf("%d %d %d %d\n",dx,dy,w,h);

  // Can disable this if you are making very big panoramas.
  // Usually this means there was an error in calculating H.
  if (false && (w > 4000 || h > 4000)) {
    printf("Can't make such big panorama :/ (%d %d)\n",w,h);
    return Image(100,100,1);
  }

  Image c(w, h, a.c);

  // Paste image a into the new image offset by dx and dy.
  for (int k = 0; k < a.c; ++k) {
    for (int j = 0; j < a.h; ++j) {
      for (int i = 0; i < a.w; ++i) {
        c(i - dx, j - dy, k) =  a(i, j, k);
      }
    }
  }

  for (int k = 0; k < b.c; ++k) {
    for (int j = topleft.y; j < botright.y; ++j) {
      for (int i = topleft.x; i < botright.x; ++i) {
        Point projected = project_point(Hba, {(double) i, (double) j});
        if (projected.y >= 0 && projected.y < b.h && projected.x >= 0 && projected.x < b.w) {
          if (j < a.h && i < a.w && i >= 0 && j >= 0) {
            if (a(i, j, k) > 0 && b(projected.x, projected.y, k) > 0) {
              float a_val = ablendcoeff * a(i, j, k);
              float temp_b;
              if (b.is_nonempty_patch(projected.x, projected.y, 1)) {
                temp_b = b.bilinear_interpolate(projected.x, projected.y, k);
              } else {
                temp_b = b(projected.x, projected.y, k);
              }
              float b_val = (1 - ablendcoeff) * temp_b;
              c.set_pixel(i - dx, j - dy, k, a_val + b_val);
            } else if (b(projected.x, projected.y, k) > 0) {
              if (b.is_nonempty_patch(projected.x, projected.y, 1)) {
                c.set_pixel(i - dx, j - dy, k, b.bilinear_interpolate(projected.x, projected.y, k));
              } else {
                c.set_pixel(i - dx, j - dy, k, b(projected.x, projected.y, k));
              }
            }
          } else {
            if (b.is_nonempty_patch(projected.x, projected.y, 1)) {
              c.set_pixel(i - dx, j - dy, k, b.bilinear_interpolate(projected.x, projected.y, k));
            } else {
              c.set_pixel(i - dx, j - dy, k, b(projected.x, projected.y, k));
            }
          }
        }
      }
    }
  }

  return trim_image(c);
}

// Create a panoramam between two images.
// const Image& a, b: images to stitch together.
// float sigma: gaussian for harris corner detector. Typical: 2
// float thresh: threshold for corner/no corner. Typical: 1-5
// int nms: window to perform nms on. Typical: 3
// float inlier_thresh: threshold for RANSAC inliers. Typical: 2-5
// int iters: number of RANSAC iterations. Typical: 1,000-50,000
// int cutoff: RANSAC inlier cutoff. Typical: 10-100
Image panorama_image(const Image& a, const Image& b, float sigma, int corner_method, float thresh, int window, int nms, float inlier_thresh, int iters, int cutoff, float acoeff) {
  // Calculate corners and descriptors
  vector<Descriptor> ad = harris_corner_detector(a, sigma, thresh, window, nms, corner_method);
  vector<Descriptor> bd = harris_corner_detector(b, sigma, thresh, window, nms, corner_method);

  // Find matches
  vector<Match> m = match_descriptors(ad, bd);

  // Run RANSAC to find the homography
  Matrix Hba = RANSAC(m, inlier_thresh, iters, cutoff);

  // Stitch the images together with the homography
  return combine_images(a, b, Hba, acoeff);
}

// Project an image onto a cylinder.
// const Image& im: image to project.
// float f: focal length used to take image (in pixels).
// returns: image projected onto cylinder, then flattened.
Image cylindrical_project(const Image& im, float f) {
  float xc = im.w / 2.0 - 1;
  float yc = im.h / 2.0 - 1;

  Image ret(im.w, im.h, im.c);
  for (int k = 0; k < ret.c; k++) {
    for (int j = 0; j < ret.h; j++) {
      for (int i = 0; i < ret.w; i++) {
        float theta = (i - xc) / f;
        float h = (j - yc) / f;

        float X_prime = sin(theta);
        float Y_prime = h;
        float Z_prime = cos(theta);

        float x_prime = f * (X_prime / Z_prime) + xc;
        float y_prime = f * (Y_prime / Z_prime) + yc;

        if (x_prime >= 0 && x_prime < im.w && y_prime >= 0 && y_prime < im.h ) {
          set_pixel(ret, i, j, k, im.bilinear_interpolate(x_prime, y_prime, k));
        }
      }
    }
  }
  return trim_image(ret);
}

// Project an image onto a cylinder.
// const Image& im: image to project.
// float f: focal length used to take image (in pixels).
// returns: image projected onto cylinder, then flattened.
Image spherical_project(const Image& im, float f) {
  float xc = im.w / 2.0 - 1;
  float yc = im.h / 2.0 - 1;

  Image ret(im.w, im.h, im.c);
  for (int k = 0; k < ret.c; k++) {
    for (int j = 0; j < ret.h; j++) {
      for (int i = 0; i < ret.w; i++) {
        float theta = (i - xc) / f;
        float sigma = (j - yc) / f;

        float X_prime = sin(theta) * cos(sigma);
        float Y_prime = sin(sigma);
        float Z_prime = cos(theta) * cos(sigma);

        float x_prime = f * (X_prime / Z_prime) + xc;
        float y_prime = f * (Y_prime / Z_prime) + yc;

        if (x_prime >= 0 && x_prime < im.w && y_prime >= 0 && y_prime < im.h ) {
          set_pixel(ret, i, j, k, im.bilinear_interpolate(x_prime, y_prime, k));
        }
      }
    }
  }
  return trim_image(ret);
}
