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

map<uint16_t, vector<Image>> scale_images(vector<Image>& im,
  uint16_t mosaicSize, uint8_t levels, bool squash) {
  map<uint16_t, vector<Image>> mapping;

  for (uint32_t level = 0; level < levels; level++) {

    int16_t internalMosaic = pow(2, level) * mosaicSize;
    vector<Image> images;

    for (unsigned i = 0; i < im.size(); i++) {
      // <parallelize>
      if (squash) {
        Image temp = bilinear_resize(im[i], internalMosaic, internalMosaic);
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

        for (uint32_t k = 0; k < im[i].c; k++) {
          for (uint32_t y = yLower; y < yUpper; y++) {
            for (uint32_t x = xLower; x < xUpper; x++) {
              square(x - xLower, y - yLower, k) = im[i](x, y, k);
            }
          }
        }
        Image temp = bilinear_resize(square, internalMosaic, internalMosaic);
        images.push_back(temp);
      }
      // </parallelize>
    }
    mapping[internalMosaic] = images;
    printf("Scaled source images to %d x %d\n", internalMosaic, internalMosaic);
  }

  return mapping;
}

Image scale_image(Image im, uint16_t mosaicSize) {

  return bilinear_resize(im, im.w / mosaicSize * mosaicSize,
    im.h / mosaicSize * mosaicSize);
}

float l1_distance(Image& a, Image& b) {
  float sum = 0.0;
  assert(a.w == b.w && a.h == b.h && a.c == b.c);
  for (uint32_t k = 0; k < a.c; k++) {
    for (uint32_t j = 0; j < a.h; j++) {
      for (uint32_t i = 0; i < a.w; i++) {
        sum += fabs(a(i, j, k) - b(i, j, k));
      }
    }
  }
  return sum;
}

float l2_distance(Image& a, Image& b) {
  float sum = 0.0;
  assert(a.w == b.w && a.h == b.h && a.c == b.c);
  for (uint32_t k = 0; k < a.c; k++) {
    for (uint32_t j = 0; j < a.h; j++) {
      for (uint32_t i = 0; i < a.w; i++) {
        sum += pow(a(i, j, k) - b(i, j, k), 2);
      }
    }
  }
  return sum;
}

vector<uint32_t> image_to_histogram(const Image& im) {
  // TODO
  vector<uint32_t> histogram;
  return histogram;
}


Image histogram_scale(Image& im, vector<uint32_t> original,
                                                    vector<uint32_t> to_match){
  // TODO
  return im;

}


uint32_t search_for_exact_match(Image& input, Image& input_dx, Image& input_dy,
  vector<Image>& source, vector<Image>& source_dx, vector<Image>& source_dy, uint8_t mode) {

  // Threshold for worse results being selected for a bit of interesing noise
  float scale = 1.7;
  // Value of derivative for select
  float derivative_scale = 0.25;

  uint32_t best_index = 0;
  float best_val = std::numeric_limits<float>::infinity();

  multimap<float, uint32_t> valueMap;

  // Sort the matches
  for (uint32_t i = 0; i < source.size(); i++) {
    float temp_val_raw = mode == 1 ? 0 : l2_distance(input, source[i]);
    float temp_val_dx = mode == 0 ? 0 : l2_distance(input_dx, source_dx[i]) * derivative_scale;
    float temp_val_dy = mode == 0 ? 0 : l2_distance(input_dy, source_dy[i]) * derivative_scale;

    float temp_val_sum = temp_val_raw + derivative_scale * (temp_val_dx + temp_val_dy);
    valueMap.insert(std::pair<float, uint32_t>(temp_val_sum, i));
    if (temp_val_sum < best_val) {
      best_index = i;
      best_val = temp_val_sum;
    }
  }

  // Pick the ones within the threshold of the best
  vector<uint32_t> candidates;
  float thresh = -1.f;
  map<float,uint32_t>::iterator it;
  for (it = valueMap.begin(); it != valueMap.end(); it++) {
    float current_val = (*it).first;
    if (thresh < 0) {
      thresh = current_val * scale;
    }
    if (current_val < thresh) {
      candidates.push_back((*it).second);
    }
  }

  return candidates[rand() % candidates.size()];
}


uint32_t search_for_fast_match(vector<uint32_t>& image_hist,
                                        vector<vector<uint32_t>> source_hist) {
  // TODO
  return 0;
}


int main(int argc, char **argv) {

  // Get input parameters
  if (argc != 9) {
    printf("USAGE: ./make-mosaic <mosaicSize> <levels> <scalePercent>");
    printf(" <squash?> <matchMethod> <sourceDir> <inputImg> <outputImg> \n");
    return 0;
  }


  // Get the base mosaic size we want
  uint16_t mosaicSize = (uint16_t) atoi(argv[1]);
  if (mosaicSize <= 0) {
    fprintf(stderr, "ERROR: mosaicSize must be > 0.");
    exit(0);
  }
  printf("Mosaic with a base tile size = %d x %d\n", mosaicSize, mosaicSize);


  // Get how many levels of different sizes we want
  int8_t levels = (int8_t) atoi(argv[2]);
  if (levels <= 0 && levels > 8) {
    fprintf(stderr, "ERROR: levels must be > 0 and <= 8.");
    exit(0);
  }
  printf("Number of larger mosaic levels = %d\n", levels);


  // Get scale factor: This determines what percent of mosaic area is made of
  // the next levels up from the current level
  double scaleFactor = atof(argv[3]);
  if (scaleFactor > 1 || scaleFactor <= 0) {
    fprintf(stderr, "ERROR: scale percent must be > 0 && < 1");
    exit(0);
  }
  printf("Level scale factor = %f\n", scaleFactor);


  // Get how we want to scale the image, 0 = crop, 1 = squash
  bool squash = (bool) atoi(argv[4]);
  printf(squash ? "Squashing source images\n" : "Cropping source images\n");


  // Get match type
  uint8_t matchMethod = (uint8_t) atoi(argv[5]);
  if (matchMethod > 3) {
    fprintf(stderr, "ERROR: match method must be 0=exact, 1=derivative, 2=exact+derivative, or 3=histogram.");
    exit(0);
  }
  printf("Match type = %s\n", matchMethod == 0 ? "Exact" : matchMethod == 1 ? "Derivative"
        : matchMethod == 2 ? "Exact + Derivative" : "Histogram");


  // Get normal source images
  vector<Image> source;
  load_images(source, argv[6]);
  printf("Loaded %lu source images\n", source.size());


  // Get input image to make into mosaic
  Image input = scale_image(load_image(string(argv[7])), mosaicSize);
  printf("Loaded %d x %d input image\n", input.w, input.h);


  // Get output file name
  string output = string(argv[8]);
  printf("Output file will be in output/%s.png\n", output.c_str());


  // Scale images into different desired sizes
  map<uint16_t, vector<Image>> mapping =
                              scale_images(source, mosaicSize, levels, squash);


  // Compute image derivatives for source images & input image
  map<uint16_t, vector<Image>> mapping_dx;
  map<uint16_t, vector<Image>> mapping_dy;
  Image input_dx;
  Image input_dy;

  if (matchMethod == 1 || matchMethod == 2) {
    Image gx_filter = make_gx_filter();
    Image gy_filter = make_gy_filter();

    for (const auto& kv : mapping) {
      vector<Image> source_dx;
      vector<Image> source_dy;
      for (uint32_t i = 0; i < kv.second.size(); i++) {
        source_dx.push_back(convolve_image(source[i], gx_filter, 1));
        source_dy.push_back(convolve_image(source[i], gy_filter, 1));
      }
      mapping_dx[kv.first] = source_dx;
      mapping_dy[kv.first] = source_dy;
      printf("Computed image derivatives for %d x %d source images\n", kv.first, kv.first);
    }

    input_dx = convolve_image(input, gx_filter, 1);
    input_dy = convolve_image(input, gy_filter, 1);
    printf("Computed image derivative for input image\n");
  }


  vector<uint16_t> scales;
  // Compute image histograms
  map<uint16_t, vector<vector<uint32_t>>> histogramMap;
  for (const auto& kv : mapping) {
    scales.push_back(kv.first);
    vector<vector<uint32_t>> histograms;
    for (uint32_t i = 0; i < kv.second.size(); i++) {
      histograms.push_back(image_to_histogram(kv.second[i]));
    }
    histogramMap[kv.first] = histograms;
    printf("Computed image histograms for %d x %d source images\n", kv.first, kv.first);
  }


  map<uint16_t, vector<uint32_t>> resultMap;
  for (auto& kv : scales) {
    printf("Processing mosaic window %d x %d\n", kv, kv);
    uint32_t processed = 0;
    uint32_t total = input.h * input.w / (kv * kv);

    for (int y = 0; y < input.h / kv; y++) {
      for (int x = 0; x < input.w / kv; x++) {
        Image input_section(kv, kv, input.c);
        Image input_section_dx(kv, kv, input.c);
        Image input_section_dy(kv, kv, input.c);
        vector<uint32_t> input_histogram = image_to_histogram(input_section);

        for (int k = 0; k < input.c; k++) {
          for (int j = 0; j < kv; j++) {
            for (int i = 0; i < kv; i++) {

              input_section(i, j, k) = input(x * kv + i, y * kv + j, k);
              if (matchMethod == 1 || matchMethod == 2) {
                input_section_dx(i, j, k) = input_dx(x * kv + i, y * kv + j, k);
                input_section_dy(i, j, k) = input_dy(x * kv + i, y * kv + j, k);
              }

            }
          }
        }

        uint32_t result = matchMethod == 3
          ? search_for_fast_match(input_histogram, histogramMap[kv])
          : search_for_exact_match(input_section, input_section_dx,
              input_section_dy, mapping[kv], mapping_dx[kv], mapping_dy[kv], matchMethod);

        processed++;
        if (!(processed % 10)) {
          printf("%d / %d processed\n", processed, total);
        }
      }
    }
  }


  // TODO: Put images back together & scale them accordingly




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


  save_png(input, "output/" + output);

  return 0;
}
