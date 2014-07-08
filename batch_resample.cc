/**
 * Batch resample the images in a directory to be isotropic (uniform
 * voxel size).
 *
 * Usage: ./batch_resample <input_dir> <output_dir> <xy_z_ratio>
 *
 */

#include <iostream>
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"

namespace fs = boost::filesystem;
typedef fs::path Path;

void ResampleImages(const Path &input_dir, const Path &output_dir, double zspacing);
void ResampleImage(const std::string &input_filename, const std::string &output_filename,
                   double zspacing);
bool EndsWith(const std::string &s, const std::string &ending);


int main (int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Usage: ./batch_resample <input_dir> <output_dir> <xy_z_ratio>" << std::endl;
    return -1;
  }

  fs::path input_dir(argv[1]);
  std::cout << input_dir << std::endl;
  fs::path output_dir(argv[2]);
  std::cout << output_dir << std::endl;
  double zspacing = atof(argv[3]);
  std::cout << zspacing << std::endl;

  if (fs::is_directory(input_dir)) {
    ResampleImages(input_dir, output_dir, zspacing);
  } else {
    std::cout << "Input dir should be a directory instead of a file. Abort." << std::endl;
  }
  return 0;
}

void ResampleImages(const Path &input_dir, const Path &output_dir, double zspacing) {
  fs::directory_iterator end_it;
  for (fs::directory_iterator it(input_dir); it != end_it; ++it) {
    std::string input_filaname = it->path().string();
    if (EndsWith(input_filaname, std::string("tif"))) {
      std::cout << input_filaname << std::endl;
      std::string output_filename = input_filaname;
      output_filename.insert(output_filename.find_last_of('.'), "-iso");
      std::cout << output_filename << std::endl;
      ResampleImage(input_filaname, output_filename, zspacing);
    }
  }
}

bool EndsWith(const std::string &s, const std::string &ending) {
  if (s.length() >= ending.length()) {
    return (0 == s.compare(s.length() - ending.length(), ending.length(), ending));
  } else {
    return false;
  }
}

void ResampleImage(const std::string &input_filename, const std::string &output_filename,
                   double zspacing) {
}
