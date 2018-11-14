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

  while ((entry = readdir(dp))) {
    if (entry != nullptr) {
        string fileName = string(entry->d_name);
        th.emplace_back(new thread([&, fileName]() {
          if (strcmp(fileName.c_str(), ".") && strcmp(fileName.c_str(), "..")
              && ends_with(fileName, string(".png"))) {
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

  printf("Loaded %d source images\n", img_total_count);
  if (!img_total_count) {
    fprintf(stderr, "ERROR: Directory contains no source images.");
    exit(0);
  }
}

void scale_images(vector<Image>& im, int16_t mosaicSize) {
  printf("%lu\n", im.size());
}

int main(int argc, char **argv) {

  // source directory, input image, mosaic size, output name
  if (argc != 5) {
    printf("USAGE: ./make-mosaic <sourceDir> <inputImg> <mosaicSize> <outputImg>\n");
    return 0;
  }

  vector<Image> im;
  load_images(im, argv[1]);

  Image srcImage = load_image(string(argv[2]));

  int16_t mosaicSize = (int16_t) atoi(argv[3]);
  if (mosaicSize <= 0) {
    fprintf(stderr, "ERROR: mosaicSize must be > 0.");
    exit(0);
  }

  string outputFile = string(argv[4]);

  scale_images(im, mosaicSize);
  // scale source images


  return 0;
}
