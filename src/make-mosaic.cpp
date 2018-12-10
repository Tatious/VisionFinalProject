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

#define WINDOW 8

bool endsWith(string const & value, string const & ending){
  if (ending.size() > value.size()) {
    return false;
  }
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

float rand_float(float min, float max) {
  return min + static_cast <float> (random()) / (RAND_MAX/(max-min));
}

void load_images(vector<Image>& im, char* indir) {
  struct dirent *entry;
  DIR *dp;

  dp = opendir(indir);
  if (dp == nullptr) {
    fprintf(stderr, "ERROR: Path does not exist or could not be read.");
    exit(0);
  }

  vector <unique_ptr<thread> > th;
  uint32_t img_total_count = 0;

  while ((entry = readdir(dp))) {
    if (entry != nullptr) {
        string fileName = string(entry->d_name);
        th.emplace_back(new thread([&, fileName]() {
          if (strcmp(fileName.c_str(), ".") && strcmp(fileName.c_str(), "..")
              && (endsWith(fileName, string(".png"))
              || endsWith(fileName, string(".jpg")))) {
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

void threadScale(bool squash, uint16_t internalMosaic, uint32_t i,
  vector<Image>& im, map<uint32_t, Image>& images) {


  if (squash) {
    Image temp = bilinear_resize(im[i], internalMosaic, internalMosaic);
    mtx.lock();
    images[i] = temp;
    mtx.unlock();
  } else {
    uint32_t sideLen = (uint32_t) min(im[i].w, im[i].h);
    Image square(sideLen, sideLen, im[i].c);

    uint32_t xLower = 0;
    uint32_t xUpper = (uint32_t) im[i].w;
    uint32_t yLower = 0;
    uint32_t yUpper = (uint32_t) im[i].h;

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
    mtx.lock();
    images[i] = temp;
    mtx.unlock();
  }
}

void threadConvolve(Image& im, Image& filter, map<uint32_t, Image>& results, uint32_t i) {
    Image res = convolve_image(im, make_gx_filter(), 1);
    mtx.lock();
    results[i] = res;
    mtx.unlock();
}

map<uint16_t, vector<Image> > scaleImages(vector<Image>& im,
  uint16_t mosaicSize, uint8_t levels, bool squash) {
  map<uint16_t, vector<Image> > mapping;

  for (uint32_t level = 0; level < levels; level++) {
    uint16_t internalMosaic = (uint16_t) pow(2, level) * mosaicSize;
    map<uint32_t, Image> images;

    vector<thread> threads;

    for (unsigned i = 0; i < im.size(); i++) {
      thread th(threadScale, squash, internalMosaic, i, std::ref(im), std::ref(images));
      threads.push_back(std::move(th));
    }

    for (auto& th : threads) th.join(); threads.clear();

    vector<Image> imagesSorted;
      for (unsigned i = 0; i < im.size(); i++) {
          imagesSorted.push_back(images[i]);
      }
    mapping[internalMosaic] = imagesSorted;
    printf("Scaled source images to %d x %d\n", internalMosaic, internalMosaic);
  }

  return mapping;
}

Image scaleImage(Image im, uint16_t mosaicSize) {

  return bilinear_resize(im, im.w / mosaicSize * mosaicSize,
    im.h / mosaicSize * mosaicSize);
}

double l1Distance(Image& a, Image& b) {
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

double l2Distance(Image& a, Image& b) {
  double sum = 0.0;
  double vals[] = {0.3, 0.6, 0.15};
  assert(a.w == b.w && a.h == b.h && a.c == b.c);
  for (uint32_t k = 0; k < a.c; k++) {

    for (uint32_t j = 0; j < a.h; j++) {
      for (uint32_t i = 0; i < a.w; i++) {
        sum += vals[k] * pow(a(i, j, k) - b(i, j, k), 2);
      }
    }
  }
  return sum;
}

Image image_to_histogram(const Image& im, const uint16_t r=8, const uint16_t g=8, const uint16_t b=8) {
  Image hist(r, g, b);

  for (uint32_t y = 0; y < im.h; y++) {
    for (uint32_t x = 0; x < im.w; x++) {
      hist((int) (im(x, y, 0) * (r - 1)), (int) (im(x, y, 1) * (g - 1)), (int) (im(x, y, 2) * (b - 1)))++;
    }
  }

  hist.feature_normalize();
  return hist;
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
    double temp_val_raw = mode == 1 ? 0 : l1Distance(input, source[i]);
    double temp_val_dx = mode == 0 || mode == 4 ? 0 : l2Distance(input_dx, source_dx[i]) * derivative_scale;
    double temp_val_dy = mode == 0 || mode == 4 ? 0 : l2Distance(input_dy, source_dy[i]) * derivative_scale;

    double temp_val_sum = temp_val_raw + derivative_scale * (temp_val_dx + temp_val_dy);

    valueMap.insert({temp_val_sum, i});
    if (temp_val_sum < best_val) {
      best_index = i;
      best_val = temp_val_sum;
    }
  }

  // Pick the ones within the threshold of the best
  vector<pair<uint32_t, double> > candidates;
  // int count = 0;
  for (auto it = valueMap.begin(); it != valueMap.end() && candidates.size() < 15; it++) {
    candidates.emplace_back((*it).second, (*it).first);
  }

  return candidates[rand() % candidates.size()];
}

pair<uint32_t, double> search_for_fast_match(Image& image_hist, vector<Image>& source_hists) {
  map<float, uint32_t, std::greater<float>> intersections;

  for (auto i = 0; i < source_hists.size(); i++) {
    Image& source_hist = source_hists[i];
    assert(source_hist.w == image_hist.w && source_hist.h == image_hist.h && source_hist.c == image_hist.c);

    double intersection = 0.0;
    for (auto c = 0; c < image_hist.c; c++) {
      for (auto y = 0; y < image_hist.h; y++) {
        for (auto x = 0; x < image_hist.w; x++) {
          intersection += fminf(image_hist(x, y, c), source_hist(x, y, c));
        }
      }
    }

    intersections.insert({intersection, i});
  }

  // Pick from the top 10% of intersections.
  vector<pair<uint32_t, double>> candidates;
  auto it = intersections.begin();
  for (auto i = 0; i < 3 && it->first != 0.0; i++, it++) {
      candidates.emplace_back(it->second, it->first);
  }

  if (candidates.empty()) {
    return {0, 0.0};
  } else {
    return candidates[rand() % candidates.size()];
  }
}

void threadMatch(uint32_t x, uint32_t y, uint32_t idx, uint32_t total,
  uint16_t scale, uint8_t matchMethod, Image& input, Image& input_dx,
  Image& input_dy, Image& outputImg, map<uint16_t, vector<Image> >& mapping,
  map<uint16_t, vector<Image> >& mapping_dx, map<uint16_t, vector<Image> >& mapping_dy,
  map<uint16_t, vector<Image>>& histogramMap, uint32_t& processed,
  vector<Image>& originalHistograms,
  vector<pair<pair<uint32_t, uint32_t>, double> >& results) {

    if (outputImg(x * scale, y * scale, 0) == -1.5) {
      Image input_section(scale, scale, input.c);
      Image input_section_dx(scale, scale, input.c);
      Image input_section_dy(scale, scale, input.c);

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

        Image input_histogram = image_to_histogram(input_section);

        Image is = matchMethod == 4
        ? bilinear_resize(input_section, WINDOW, WINDOW)
        : input_section;
        scale = matchMethod == 4 ? WINDOW : scale;
        pair<uint32_t, double> result = matchMethod == 3
        ? search_for_fast_match(input_histogram, histogramMap[scale])
        : search_for_exact_match(is,
                input_section_dx,
            input_section_dy,
            mapping[scale],
            mapping_dx[scale],
            mapping_dy[scale],
            matchMethod);

      pair<uint32_t, uint32_t> p1(idx, result.first);
      pair<pair<uint32_t, uint32_t>, double> p2(p1, result.second);

      mtx.lock();
      originalHistograms.push_back(std::move(input_histogram));
      results.push_back(p2);
      mtx.unlock();
    }

    mtx.lock();
    processed++;
    if (!(processed % 128)) {
      printf("%d / %d processed\n", processed, total);
    }
    mtx.unlock();

}

Image meanScale(Image& im1, Image& im2, float meanFactor) {
    Image im3(im1.w, im1.h, im1.c);
    for (uint32_t k = 0; k < im1.c; k++) {
        double m1 = 0.0;
        double m2 = 0.0;
        uint32_t count = 0;

        for (uint32_t j = 0; j < im1.h; j++) {
            for (uint32_t i = 0; i < im1.w; i++) {
                m1 += im1(i, j, k);
                m2 += im2(i, j, k);
                count++;
            }
        }
        m1 /= count;
        m2 /= count;
        double m3 = (m1 - m2) * meanFactor;
        for (uint32_t j = 0; j < im1.h; j++) {
            for (uint32_t i = 0; i < im1.w; i++) {
                im3(i, j, k) = im1(i, j, k) - m3;
            }
        }
    }

    return im3;
}

void threadCombine(Image& outputImg, Image& toUse, Image& original, Image& origHist,
  Image& matchHist, uint16_t scale, float meanFactor, uint32_t inputC, uint32_t xStart,
  uint32_t yStart) {

    Image a(toUse.w, toUse.h, toUse.c);
    for (int k = 0; k < inputC; k++) {
        for (int j = 0; j < scale; j++) {
            for (int i = 0; i < scale; i++) {
                a(i, j, k) = original(xStart + i, yStart + j, k);
            }
        }
    }
    Image b = meanScale(toUse, a, meanFactor);

    for (int k = 0; k < inputC; k++) {
      for (int j = 0; j < scale; j++) {
        for (int i = 0; i < scale; i++) {
          set_pixel(outputImg, xStart + i, yStart + j, k, b(i, j, k));
        }
      }
    }
}






int main(int argc, char **argv) {

  vector<thread> threads;

  // Get input parameters
  if (argc != 10) {
    printf("USAGE: ./make-mosaic <mosaicSize> <levels> <scalePercent> <meanAmount>");
    printf(" <squash?> <matchMethod> <sourceDir> <inputImg> <outputImg> \n");
    return 0;
  }


  // Get the base mosaic size we want
  auto mosaicSize = (uint16_t) atoi(argv[1]);
  if (mosaicSize <= 0) {
    fprintf(stderr, "ERROR: mosaicSize must be > 0.");
    return 1;
  }
  printf("Mosaic with a base tile size = %d x %d\n", mosaicSize, mosaicSize);


  // Get how many levels of different sizes we want
  auto levels = (uint8_t) atoi(argv[2]);
  if (levels > 8) {
    fprintf(stderr, "ERROR: levels must be <= 8.");
    return 1;
  }
  printf("Number of larger mosaic levels = %d\n", levels);


  // Get scale factor: This determines what percent of mosaic area is made of
  // the next levels up from the current level
  double scaleFactor = atof(argv[3]);
  if (scaleFactor > 1 || scaleFactor <= 0) {
    fprintf(stderr, "ERROR: scale percent must be > 0 && < 1");
    return 1;
  }
  printf("Level scale factor = %f\n", scaleFactor);


  float meanFactor = atof(argv[4]);
  if (meanFactor <= 0) {
    fprintf(stderr, "ERROR: mean factor must be > 0");
    return 1;
  }


  // Get how we want to scale the image, 0 = crop, 1 = squash
  bool squash = (bool) atoi(argv[5]);
  printf(squash ? "Squashing source images\n" : "Cropping source images\n");


  // Get match type
  auto matchMethod = (uint8_t) atoi(argv[6]);
  if (matchMethod > 5) {
    fprintf(stderr, "ERROR: match method must be 0=exact, 1=derivative, 2=exact+derivative, 3=histogram, or 4=fast exact");
    exit(0);
  }
  printf("Match type = %s\n", matchMethod == 0 ? "Exact" : matchMethod == 1 ? "Derivative"
        : matchMethod == 2 ? "Exact + Derivative" : matchMethod == 3 ? "Histogram"
      : "Fast Exact");


  // Get normal source images
  vector<Image> source;
  load_images(source, argv[7]);
  printf("Loaded %lu source images\n", source.size());


  // Get input image to make into mosaic
  vector<Image> inputImages;
  if (true) {
    inputImages.push_back(load_image(string(argv[8])));
  } else {
    load_images(inputImages, argv[8]);
  }
  for(int32_t i = 0; i < inputImages.size(); i++) {
    inputImages[i] = bilinear_resize(inputImages[i], inputImages[i].w / mosaicSize * mosaicSize, inputImages[i].h / mosaicSize * mosaicSize);
  }
  printf("Loaded %lu input images\n", inputImages.size());


  // Get output file name
  string output = string(argv[9]);
  printf("Output file will be in output/%s.png\n", output.c_str());


  // Scale images into different desired sizes
  map<uint16_t, vector<Image> > mapping = scaleImages(source, mosaicSize, levels, squash);

  map<uint16_t, vector<Image> > matchMapping = matchMethod == 4
  ? scaleImages(source, WINDOW, 1, squash)
  : mapping;

  // Compute image derivatives for source images & input image
  map<uint16_t, vector<Image> > mappingDx;
  map<uint16_t, vector<Image> > mappingDy;
  vector<Image> input_dx_images;
  vector<Image> input_dy_images;

  if (matchMethod == 1 || matchMethod == 2) {
    Image gxFilter = make_gx_filter();
    Image gyFilter = make_gy_filter();

    map<uint32_t, Image> imDx;
    map<uint32_t, Image> imDy;
    vector<thread> threads;
    for (uint32_t i = 0; i < source.size(); i++) {
        thread th(threadConvolve, std::ref(source[i]), std::ref(gxFilter), std::ref(imDx), i);
        threads.push_back(std::move(th));
        th = thread(threadConvolve, std::ref(source[i]), std::ref(gyFilter), std::ref(imDy), i);
        threads.push_back(std::move(th));
    }
    for (auto& th : threads) th.join(); threads.clear();

    vector<Image> sourceDx;
    vector<Image> sourceDy;
    for (uint32_t i = 0; i < source.size(); i++) {
        sourceDx.push_back(imDx[i]);
        sourceDy.push_back(imDy[i]);
    }
    mappingDx = scaleImages(sourceDx, mosaicSize, levels, squash);
    mappingDy = scaleImages(sourceDy, mosaicSize, levels, squash);

    for (auto& input : inputImages) {
      input_dx_images.push_back(convolve_image(input, gxFilter, 1));
      input_dy_images.push_back(convolve_image(input, gyFilter, 1));
    }
    printf("Computed image derivative for input image\n");
  } else {
    for (const auto& kv : mapping) {
      vector<Image> source_dx;
      vector<Image> source_dy;
      mappingDx[kv.first] = source_dx;
      mappingDy[kv.first] = source_dy;
    }
  }




  vector<uint16_t> scales;
  // Compute image histograms
  map<uint16_t, vector<Image>> histogramMap;
  for (const auto& kv : mapping) {
    scales.push_back(kv.first);
    vector<Image> histograms;
    for (const auto &i : kv.second) {
      histograms.push_back(image_to_histogram(i));
    }
    histogramMap.insert({kv.first, std::move(histograms)});
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

  // Do this backwards to avoid recomputation
  for (int32_t i = 0; i < inputImages.size(); i++) {
    Image& input = inputImages[i];
    Image& input_dx = input_dx_images.empty() ? input : input_dx_images[i];
    Image& input_dy = input_dx_images.empty() ? input : input_dy_images[i];
    // Initialize as negative so we know what's been filled in
    Image outputImg(input.w, input.h, input.c);
    for (int k = 0; k < outputImg.c; k++) {
      for (int y = 0; y < outputImg.h; y++) {
        for (int x = 0; x < outputImg.w; x++) {
          outputImg(x, y, k) = -1.5f;
        }
      }
    }

    double imageSize = input.w * input.h;
    for (auto i = (int32_t) (scales.size() - 1); i >= 0; i--) {
      vector<pair<pair<uint32_t, uint32_t>, double> > results; // location, match index, score
      vector<Image> originalHistograms;
      uint16_t scale = scales[i];
      printf("Processing mosaic window %d x %d...\n", scale, scale);
      uint32_t processed = 0;
      auto total = (uint32_t) ((input.h * input.w) / (scale * scale));

      uint32_t idx = 0;
      for (uint32_t y = 0; y < input.h / scale; y++) {
        for (uint32_t x = 0; x < input.w / scale; x++) {

          thread th(threadMatch, x, y, idx, total, scale, matchMethod,
            std::ref(input), std::ref(input_dx), std::ref(input_dy),
            std::ref(outputImg), std::ref(matchMapping), std::ref(mappingDx),
            std::ref(mappingDy), std::ref(histogramMap), std::ref(processed),
            std::ref(originalHistograms), std::ref(results));
          threads.push_back(std::move(th));

          idx++;
        }
      }

      for (auto& th : threads)th.join();threads.clear();


      printf("Using %.2f%% of these bad bois\n", percentages[i] * 100);

      double segmentSize = scale * scale;
      int cnt = (int) ceil((imageSize * percentages[i]) / segmentSize);
      printf("ADDING %d OF THESE\n", cnt);

      // Sort the results
      std::sort(results.begin(), results.end(), [](const std::pair<pair<uint32_t, uint32_t>, double> &left,
        const std::pair<pair<uint32_t, uint32_t>, double> &right) {
        return left.second < right.second;
      });

      int xBlocks = input.w / scale;
      int count = 0;
      idx = 0;

      threads.clear();
      while (count < cnt && idx < results.size()) {
        // These are essentially 0,0: we move on in scale X scale chunks

        int xBlk = results[idx].first.first % xBlocks;
        int yBlk = results[idx].first.first / xBlocks;

        int xStart = xBlk * scale;
        int yStart = yBlk * scale;

        if (outputImg(xStart, yStart, 0) == -1.5) {
          Image& toUse = mapping[scale][results[idx].first.second];
          Image& origHist = originalHistograms[idx];
          Image& matchHist = histogramMap[scale][results[idx].first.second];

          thread th(threadCombine, std::ref(outputImg), std::ref(toUse), std::ref(input), std::ref(origHist),
            std::ref(matchHist), scale, meanFactor, input.c, xStart, yStart);
          threads.push_back(std::move(th));

          count++;
        }
        idx++;
      }

      for (auto& th : threads)th.join();threads.clear();

    }

    printf("%d/%lu\n", i+1, inputImages.size());
    outputImg.clamp();
    //save_png(outputImg, string(argv[8]) + "/result" + to_string(i));
    save_png(outputImg, "output/" + output);
  }

  return 0;
}
