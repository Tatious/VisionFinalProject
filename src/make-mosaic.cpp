#include "image.h"
#include "utils.h"
#include "matrix.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include <string>
#include <thread>
#include <map>
#include <mutex>

using namespace std;

std::mutex mtx;

bool ends_with(string const & value, string const & ending){
  if (ending.size() > value.size()) {
    return false;
  }
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

void load_images(vector<Image>& im, char* indir) {
  struct dirent *entry;
  DIR *dp;

  dp = opendir(indir);
  if (dp == NULL) {
    fprintf(stderr, "ERROR: Path does not exist or could not be read.");
    exit(0);
  }

  vector <unique_ptr<thread>> th;
  uint32_t img_total_count = 0;

  // TODO: Make parallelization faster
  while ((entry = readdir(dp))) {
    if (entry != nullptr) {
        string fileName = string(entry->d_name);
        th.emplace_back(new thread([&, fileName]() {
          if (strcmp(fileName.c_str(), ".") && strcmp(fileName.c_str(), "..")
              && (ends_with(fileName, string(".png"))
              || ends_with(fileName, string(".jpg")))) {
            string file = string(indir) + "/" + fileName;
            Image i = load_image(file);
            mtx.lock();
            img_total_count++;
            im.push_back(i);
            mtx.unlock();
          }
        }));
    }
  }

  closedir(dp);

  for (auto&e1:th)e1->join();th.clear();

  printf("Loaded %d raw source images\n", img_total_count);
  if (!img_total_count) {
    fprintf(stderr, "ERROR: Directory contains no source images.");
    exit(0);
  }
}

void scale_images(vector<Image>& im, int16_t mosaicSize, bool squash) {
  // TODO: Parallelize
  for (unsigned i=0; i < im.size(); i++) {
    if (squash) {
      im[i] = bilinear_resize(im[i], mosaicSize, mosaicSize);
    } else {
      uint32_t sideLen = min(im[i].w, im[i].h);
      Image square(sideLen, sideLen, im[i].c);

      uint32_t xLower = 0;
      uint32_t xUpper = im[i].w;
      uint32_t yLower = 0;
      uint32_t yUpper = im[i].h;

      if (im[i].w != sideLen) {
        xLower += (uint32_t) floor((xUpper - sideLen) / 2.f);
        xUpper -= (uint32_t) floor((xUpper - sideLen) / 2.f);
      }
      if (im[i].h != sideLen) {
        yLower += (uint32_t) floor((yUpper - sideLen) / 2.f);
        yUpper -= (uint32_t) floor((yUpper - sideLen) / 2.f);
      }

      for (int k = 0; k < im[i].c; k++) {
        for (int y = yLower; y < yUpper; y++) {
          for (int x = xLower; x < xUpper; x++) {
            square(x, y, k) = im[i](x, y, k);
          }
        }
      }
      im[i] = bilinear_resize(square, mosaicSize, mosaicSize);
    }
  }
  printf("Scaled source images to %d x %d\n", mosaicSize, mosaicSize);
}

Image scale_image(Image im, int16_t mosaicSize) {

  return bilinear_resize(im, im.w / mosaicSize * mosaicSize,
    im.h / mosaicSize * mosaicSize);
}

int main(int argc, char **argv) {
  if (argc != 6) {
    printf("USAGE: ./make-mosaic <sourceDir> <inputImg> <mosaicSize>");
    printf(" <squash?> <outputImg>\n");
    return 0;
  }
  string outputFile = string(argv[5]);

  vector<Image> im;
  load_images(im, argv[1]);

  int16_t mosaicSize = (int16_t) atoi(argv[3]);
  if (mosaicSize <= 0) {
    fprintf(stderr, "ERROR: mosaicSize must be > 0.");
    exit(0);
  }
  bool squash = (bool) atoi(argv[4]);

  Image inputImage = scale_image(load_image(string(argv[2])), mosaicSize);
  printf("Loaded %d x %d input image\n", inputImage.w, inputImage.h);

  scale_images(im, mosaicSize, squash);

  // Split srcImage up into subimages & search best match foreach of them
  // Compare:
  //    raw pixel values & derivative of image + scaled in RGB (binary search)

  return 0;
}
