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
#include <limits>

using namespace std;

std::mutex mtx;

bool ends_with(string const & value, string const & ending){
  if (ending.size() > value.size()) {
    return false;
  }
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

float rand_float(float min, float max) {
  return min + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(max-min)));
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

map<uint8_t, vector<Image>> scale_images(vector<Image>& im,
  int16_t mosaicSize, uint8_t levels, bool squash) {
  map<uint8_t, vector<Image>> mapping;

  for (int level = 0; level < levels; level++) {

    int16_t internalMosaic = pow(2, level) * mosaicSize;
    vector<Image> images;

    for (unsigned i = 0; i < im.size(); i++) {
      // <parallelize>
      if (squash) {
        Image temp = bilinear_resize(im[i], internalMosaic, internalMosaic);
        rgb_to_hsv(temp);
        images.push_back(temp);
      } else {
        uint32_t sideLen = min(im[i].w, im[i].h);
        Image square(sideLen, sideLen, im[i].c);

        uint32_t xLower = 0;
        uint32_t xUpper = im[i].w;
        uint32_t yLower = 0;
        uint32_t yUpper = im[i].h;

        if (im[i].w != sideLen) {
          xLower += (uint32_t) floor((xUpper - sideLen) / 2.f);
          xUpper -= (uint32_t) ceil((xUpper - sideLen) / 2.f);
        }
        if (im[i].h != sideLen) {
          yLower += (uint32_t) floor((yUpper - sideLen) / 2.f);
          yUpper -= (uint32_t) ceil((yUpper - sideLen) / 2.f);
        }

        for (int k = 0; k < im[i].c; k++) {
          for (int y = yLower; y < yUpper; y++) {
            for (int x = xLower; x < xUpper; x++) {
              square(x - xLower, y - yLower, k) = im[i](x, y, k);
            }
          }
        }
        Image temp = bilinear_resize(square, internalMosaic, internalMosaic);
        rgb_to_hsv(temp);
        images.push_back(temp);
      }
      // </parallelize>
    }
    mapping[internalMosaic] = images;
    printf("Scaled source images to %d x %d\n", internalMosaic, internalMosaic);
  }

  return mapping;
}

Image scale_image(Image im, int16_t mosaicSize) {

  return bilinear_resize(im, im.w / mosaicSize * mosaicSize,
    im.h / mosaicSize * mosaicSize);
}

float l1_distance(Image& a, Image& b) {
  float sum = 0.0;
  assert(a.w == b.w && a.h == b.h && a.c == b.c);
  for (int k = 0; k < a.c; k++) {
    for (int j = 0; j < a.h; j++) {
      for (int i = 0; i < a.w; i++) {
        sum += fabs(a(i, j, k) - b(i, j, k));
      }
    }
  }
  return sum;
}

float l2_distance(Image& a, Image& b) {
  float sum = 0.0;
  assert(a.w == b.w && a.h == b.h && a.c == b.c);
  for (int k = 0; k < a.c; k++) {
    for (int j = 0; j < a.h; j++) {
      for (int i = 0; i < a.w; i++) {
        sum += pow(a(i, j, k) - b(i, j, k), 2);
      }
    }
  }
  return sum;
}

Image search_for_match(Image& input, Image& input_dx, Image& input_dy,
  vector<Image>& source, vector<Image>& source_dx, vector<Image>& source_dy) {

  float scale = 1.7;
  float derivative_scale = 0.25;

  uint32_t best_index = 0;
  float best_val = std::numeric_limits<float>::infinity();

  multimap<float, int> valueMap;

  // Sort the matches
  for (int i = 0; i < source.size(); i++) {
    float temp_val_raw = l2_distance(input, source[i]);
    float temp_val_dx = l2_distance(input_dx, source_dx[i]) * derivative_scale;
    float temp_val_dy = l2_distance(input_dy, source_dy[i]) * derivative_scale;

    float temp_val_sum = temp_val_raw + temp_val_dx + temp_val_dy;
    valueMap.insert(std::pair<float, int>(temp_val_sum, i));
    if (temp_val_sum < best_val) {
      best_index = i;
      best_val = temp_val_sum;
    }
  }

  // Pick the ones within the threshold of the best
  vector<int> candidates;
  float thresh = -1.f;
  map<float,int> :: iterator it;
  for (it=valueMap.begin() ; it!=valueMap.end() ; it++) {
    float current_val = (*it).first;
    if (thresh < 0) {
      thresh = current_val * scale;
    }
    if (current_val < thresh) {
      candidates.push_back((*it).second);
    }
  }
  Image best_match = source[candidates[rand() % candidates.size()]];

  // Scale the Hue
  float diff = 0.0;
  float sat = 0.0;
  for (int j = 0; j < input.h; j++) {
    for (int i = 0; i < input.w; i++) {
      diff += best_match(i, j, 0) == 0 ? 0 : input(i, j, 0) / best_match(i, j, 0);
      sat += input(i, j, 1);
    }
  }

  scale_image(best_match, 0, diff);

  return best_match;
}

int main(int argc, char **argv) {

  if (argc != 7) {
    printf("USAGE: ./make-mosaic <sourceDir> <inputImg> <outputImg> <mosaicSize>");
    printf(" <levels> <squash?> \n");
    return 0;
  }
  string output = string(argv[3]);

  vector<Image> source;
  load_images(source, argv[1]);

  int16_t mosaicSize = (int16_t) atoi(argv[4]);
  if (mosaicSize <= 0) {
    fprintf(stderr, "ERROR: mosaicSize must be > 0.");
    exit(0);
  }

  int8_t levels = (int8_t) atoi(argv[5]);
  if (levels <= 0 && levels > 8) {
    fprintf(stderr, "ERROR: levels must be > 0 and <= 8.");
    exit(0);
  }
  bool squash = (bool) atoi(argv[6]);

  Image input = scale_image(load_image(string(argv[2])), mosaicSize);
  rgb_to_hsv(input);
  printf("Loaded %d x %d input image\n", input.w, input.h);

  map<uint8_t, vector<Image>> mapping =
                              scale_images(source, mosaicSize, levels, squash);


  map<uint8_t, vector<Image>> mapping_dx;
  map<uint8_t, vector<Image>> mapping_dy;
  Image gx_filter = make_gx_filter();
  Image gy_filter = make_gy_filter();

  for (const auto& kv : mapping) {
    vector<Image> source_dx;
    vector<Image> source_dy;
    for (int i = 0; i < kv.second.size(); i++) {
      source_dx.push_back(convolve_image(source[i], gx_filter, 1));
      source_dy.push_back(convolve_image(source[i], gy_filter, 1));
    }
    mapping_dx[kv.first] = source_dx;
    mapping_dy[kv.first] = source_dy;
  }



  Image input_dx = convolve_image(input, gx_filter, 1);
  Image input_dy = convolve_image(input, gy_filter, 1);
  printf("Computed image derivatives\n");


  uint32_t processed = 0;
  uint32_t total = input.h * input.w / (mosaicSize * mosaicSize);
  //vector <unique_ptr<thread>> th;

  /*
  for (int y = 0; y < input.h / mosaicSize; y++) {
    for (int x = 0; x < input.w / mosaicSize; x++) {
      // <parallelize>

      //th.emplace_back(new thread([&, x, y]() {
      Image input_section(mosaicSize, mosaicSize, input.c);
      Image input_section_dx(mosaicSize, mosaicSize, input.c);
      Image input_section_dy(mosaicSize, mosaicSize, input.c);
      for (int k = 0; k < input.c; k++) {
        for (int j = 0; j < mosaicSize; j++) {
          for (int i = 0; i < mosaicSize; i++) {
            input_section(i, j, k) = input(x * mosaicSize + i, y * mosaicSize + j, k);
            input_section_dx(i, j, k) = input_dx(x * mosaicSize + i, y * mosaicSize + j, k);
            input_section_dy(i, j, k) = input_dy(x * mosaicSize + i, y * mosaicSize + j, k);
          }
        }
      }

      Image result = search_for_match(input_section, input_section_dx,
        input_section_dy, source, source_dx, source_dy);

      for (int k = 0; k < input.c; k++) {
        for (int j = 0; j < mosaicSize; j++) {
          for (int i = 0; i < mosaicSize; i++) {
            input(x * mosaicSize + i, y * mosaicSize + j, k) = result(i, j, k);
          }
        }
      }
      //mtx.lock();
      processed++;
      if (!(processed % 10)) {
        printf("%d / %d processed\n", processed, total);
      }

      //mtx.unlock();
      //}));
      // </parallelize>
    }
  }
  */
  
  //for (auto&e1:th)e1->join();th.clear();
  hsv_to_rgb(input);
  save_png(input, "output/" + output);

  return 0;
}
