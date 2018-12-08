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

  vector <unique_ptr<thread> > th;
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

map<uint16_t, vector<Image> > scale_images(vector<Image>& im,
  uint16_t mosaicSize, uint8_t levels, bool squash) {
  map<uint16_t, vector<Image> > mapping;

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

double l1_distance(Image& a, Image& b) {
  double sum = 0.0;
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

double l2_distance(Image& a, Image& b) {
  double sum = 0.0;
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


pair<uint32_t, double> search_for_exact_match(Image& input, Image& input_dx, Image& input_dy,
  vector<Image>& source, vector<Image>& source_dx, vector<Image>& source_dy, uint8_t mode) {

  // Threshold for worse results being selected for a bit of interesing noise
  double scale = 1.7;
  // Value of derivative for select
  double derivative_scale = 0.25;

  uint32_t best_index = 0;
  double best_val = std::numeric_limits<float>::infinity();

  multimap<double, uint32_t> valueMap;
  // Sort the matches
  for (uint32_t i = 0; i < source.size(); i++) {
    double temp_val_raw = mode == 1 ? 0 : l2_distance(input, source[i]);
    double temp_val_dx = mode == 0 ? 0 : l2_distance(input_dx, source_dx[i]) * derivative_scale;
    double temp_val_dy = mode == 0 ? 0 : l2_distance(input_dy, source_dy[i]) * derivative_scale;

    double temp_val_sum = temp_val_raw + derivative_scale * (temp_val_dx + temp_val_dy);
    valueMap.insert(pair<double, uint32_t>(temp_val_sum, i));
    if (temp_val_sum < best_val) {
      best_index = i;
      best_val = temp_val_sum;
    }
  }

  // Pick the ones within the threshold of the best
  vector<pair<uint32_t, double> > candidates;
  double thresh = -1.f;
  map<double,uint32_t>::iterator it;
  for (it = valueMap.begin(); it != valueMap.end(); it++) {
    double current_val = (*it).first;
    if (thresh < 0) {
      thresh = current_val * scale;
    }
    if (current_val < thresh) {
      candidates.push_back(pair<uint32_t, double>((*it).second, (*it).first));
    }
  }

  return candidates[rand() % candidates.size()];
}


pair<uint32_t, double> search_for_fast_match(vector<uint32_t>& image_hist,
                                        vector<vector<uint32_t> > source_hist) {
  // TODO
  return pair<uint32_t, double> (0, 0.0);
}

void threadMatch(uint32_t x, uint32_t y, uint32_t idx, uint32_t total,
  uint16_t scale, uint8_t matchMethod, Image& input, Image& input_dx,
  Image& input_dy, Image& outputImg, map<uint16_t, vector<Image> >& mapping,
  map<uint16_t, vector<Image> >& mapping_dx, map<uint16_t, vector<Image> >& mapping_dy,
  map<uint16_t, vector<vector<uint32_t> > >& histogramMap, uint32_t& processed,
  vector<vector<uint32_t> >& originalHistograms,
  vector<pair<pair<uint32_t, uint32_t>, double> >& results) {

    if (outputImg(x * scale, y * scale, 0) == -1.5) {
      Image input_section(scale, scale, input.c);
      Image input_section_dx(scale, scale, input.c);
      Image input_section_dy(scale, scale, input.c);
      vector<uint32_t> input_histogram = image_to_histogram(input_section);

      mtx.lock();
      originalHistograms.push_back(input_histogram);
      mtx.unlock();

      for (int k = 0; k < input.c; k++) {
        for (int j = 0; j < scale; j++) {
          for (int i = 0; i < scale; i++) {
            input_section(i, j, k) = input(x * scale + i, y * scale + j, k);
            if (matchMethod == 1 || matchMethod == 2) {
              input_section_dx(i, j, k) = input_dx(x * scale + i, y * scale + j, k);
              input_section_dy(i, j, k) = input_dy(x * scale + i, y * scale + j, k);
            }
          }
        }
      }

      pair<uint32_t, double> result = matchMethod == 3
        ? search_for_fast_match(input_histogram, histogramMap[scale])
        : search_for_exact_match(input_section, input_section_dx,
            input_section_dy, mapping[scale], mapping_dx[scale], mapping_dy[scale], matchMethod);

      pair<uint32_t, uint32_t> p1(idx, result.first);
      pair<pair<uint32_t, uint32_t>, double> p2(p1, result.second);

      mtx.lock();
      results.push_back(p2);
      mtx.unlock();
    }

    mtx.lock();
    processed++;
    if (!(processed % 10)) {
      printf("%d / %d processed\n", processed, total);
    }
    mtx.unlock();

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
  map<uint16_t, vector<Image> > mapping =
                              scale_images(source, mosaicSize, levels, squash);


  // Compute image derivatives for source images & input image
  map<uint16_t, vector<Image> > mapping_dx;
  map<uint16_t, vector<Image> > mapping_dy;
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
  map<uint16_t, vector<vector<uint32_t> > > histogramMap;
  for (const auto& kv : mapping) {
    scales.push_back(kv.first);
    vector<vector<uint32_t> > histograms;
    for (uint32_t i = 0; i < kv.second.size(); i++) {
      histograms.push_back(image_to_histogram(kv.second[i]));
    }
    histogramMap[kv.first] = histograms;
    printf("Computed image histograms for %d x %d source images\n", kv.first, kv.first);
  }


  // Compute the percentages for each scale size to use
  vector<double> percentages;
  double percentLeft = 1.0;
  for (uint8_t i = 0; i < levels - 1; i++) {
    percentages.push_back(percentLeft * (1 - scaleFactor));
    percentLeft *= scaleFactor;
  }
  percentages.push_back(percentLeft);


  // Initialize as negative so we know what's been filled in
  Image outputImg(input.w, input.h, input.c);
  for (int k = 0; k < outputImg.c; k++) {
    for (int y = 0; y < outputImg.h; y++) {
      for (int x = 0; x < outputImg.w; x++) {
        outputImg(x, y, k) = -1.5;
      }
    }
  }

  // Do this backwards to avoid recomputation
  double imageSize = input.w * input.h;
  for (int i = scales.size() - 1; i >= 0; i--) {
    vector<pair<pair<uint32_t, uint32_t>, double> > results; // location, match index, score
    vector<vector<uint32_t> > originalHistograms;
    uint16_t scale = scales[i];
    printf("Processing mosaic window %d x %d...\n", scale, scale);
    uint32_t processed = 0;
    uint32_t total = input.h * input.w / (scale * scale);

    std::thread t[(input.h / scale) * (input.w / scale)];
    uint32_t idx = 0;
    for (uint32_t y = 0; y < input.h / scale; y++) {
      for (uint32_t x = 0; x < input.w / scale; x++) {

        t[idx] = thread(threadMatch, x, y, idx, total, scale, matchMethod,
          std::ref(input), std::ref(input_dx), std::ref(input_dy),
          std::ref(outputImg), std::ref(mapping), std::ref(mapping_dx),
          std::ref(mapping_dy), std::ref(histogramMap), std::ref(processed),
          std::ref(originalHistograms), std::ref(results));

        idx++;
      }
    }

    for (int thr = 0; thr < (input.h / scale) * (input.w / scale); thr++) {
      t[thr].join();
    }

    printf("Using %.2f%% of these bad bois\n", percentages[i] * 100);

    double segmentSize = scale * scale;
    int cnt = (int) ceil((imageSize * percentages[i]) / segmentSize);
    printf("ADDING %d OF THESE\n", cnt);

    // Sort the results
    std::sort(results.begin(), results.end(), [](const std::pair<pair<uint32_t, uint32_t>, double> &left,
      const std::pair<pair<uint32_t, uint32_t>, double> &right) {
      return left.second < right.second;
    });

    int x_blocks = input.w / scale;
    int y_blocks = input.h / scale;
    int count = 0;
    idx = 0;
    while (count < cnt && idx < results.size()) {
      // These are essentially 0,0: we move on in scale X scale chunks

      int x_blk = results[idx].first.first % x_blocks;
      int y_blk = results[idx].first.first / x_blocks;

      int x_start = x_blk * scale;
      int y_start = y_blk * scale;

      if (outputImg(x_start, y_start, 0) == -1.5) {
        Image toUse = mapping[scale][results[idx].first.second];
        vector<uint32_t> origHist = originalHistograms[idx];
        vector<uint32_t> matchHist = histogramMap[scale][results[idx].first.second];

        toUse = histogram_scale(toUse, origHist, matchHist);

        for (int k = 0; k < input.c; k++) {
          for (int j = 0; j < scale; j++) {
            for (int i = 0; i < scale; i++) {
              set_pixel(outputImg, x_start + i, y_start + j, k, get_pixel(toUse, i, j, k));
            }
          }
        }
        count++;
      }
      idx++;
    }
  }


  save_png(outputImg, "output/" + output);

  return 0;
}
