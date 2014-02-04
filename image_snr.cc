#include <iostream>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <boost/filesystem.hpp>
#include "multisnake.h"


/*
 * Compute the image SNR using Otsu's binary segmentation.
 */

std::string GetSuffix(const std::string &s) {
  return s.substr(s.size()-3);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: ./snr <image_path>" << std::endl;
    return EXIT_FAILURE;
  }

  namespace fs = boost::filesystem;
  fs::path image_dir(argv[1]);

  try {
    if (fs::exists(image_dir)) {
      soax::Multisnake multisnake;
      // fs::directory_iterator end_it;
      // for (fs::directory_iterator it(image_dir); it != end_it; ++it) {
      //   multisnake.LoadImage(it->path().string());
      //   double snr = multisnake.ComputeImageSNR();
      //   std::cout << it->path().filename() << " SNR: " << snr << std::endl;
      // }
      typedef std::vector<fs::path> Paths;
      Paths image_paths;
      std::copy(fs::directory_iterator(image_dir), fs::directory_iterator(),
                back_inserter(image_paths));
      std::sort(image_paths.begin(), image_paths.end());
      for (Paths::const_iterator it(image_paths.begin());
           it != image_paths.end(); ++it) {
        if (GetSuffix(it->string()) != "mha") continue;
        multisnake.LoadImage(it->string());
        double snr = multisnake.ComputeImageSNR();
        // if (snr < 4.5 && snr > 3.5)
        std::cout << it->filename() << " SNR: " << snr << std::endl;
      }
    } else {
      std::cout << image_dir << " does not exist." << std::endl;
    }
  } catch (const fs::filesystem_error &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
