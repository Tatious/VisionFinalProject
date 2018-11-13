#include "image.h"
#include "utils.h"
#include "matrix.h"

#include <string>
#include <thread>
#include <map>
#include <mutex>

using namespace std;


int main(int argc, char **argv)
  {

  // source directory, input image, mosaic size, output name
  if(argc != 5)
    {
    printf("USAGE: ./make-mosaic <sourceDir> <inputImg> <mosaicSize> <outputImg>\n");
    return 0;
    }

  // argv[1] : source directory
  // argv[2] : input image
  // argv[3] : mosaic size
  // argv[4] : outputName

  return 0;
  }
