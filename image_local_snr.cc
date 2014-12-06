/*
 * Batch compute local image SNRs using ground truth SOACs.
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <boost/filesystem.hpp>
#include "global.h"
#include "utility.h"
#include "multisnake.h"


std::string GetSuffix(const std::string &s) {
  return s.substr(s.size()-3);
}



int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Usage: ./snr <image_dir> <snake_dir> <outfile_path>" << std::endl;
    return EXIT_FAILURE;
  }

  namespace fs = boost::filesystem;
  fs::path image_dir(argv[1]);

  try {
    if (fs::exists(image_dir)) {
      soax::Multisnake multisnake;

      std::ofstream outfile(argv[3]);
      if (!outfile) {
        std::cerr << "could not open file: " << argv[3] << std::endl;
        return EXIT_FAILURE;
      }

      typedef std::vector<fs::path> Paths;
      Paths image_paths;
      std::copy(fs::directory_iterator(image_dir), fs::directory_iterator(),
                back_inserter(image_paths));
      std::sort(image_paths.begin(), image_paths.end());
      for (Paths::const_iterator it(image_paths.begin());
           it != image_paths.end(); ++it) {
        // if (fs::is_regular_file(image_dir)) {
        if (GetSuffix(it->string()) == "mha" ||
            GetSuffix(it->string()) == "tif") {
        // if (GetSuffix(image_dir.string()) == "mha" ||
        //     GetSuffix(image_dir.string()) == "tif") {
          multisnake.LoadImage(it->string());
          std::cout << it->string() << std::endl;
          // multisnake.LoadImage(image_dir.string());
          std::string gt_path = std::string(argv[2]) + it->filename().string();
          gt_path.replace(gt_path.size() - 3, 3, "txt");
          std::cout << gt_path << std::endl;
          multisnake.LoadGroundTruthSnakes(gt_path);
          std::cout << multisnake.GetNumberOfComparingSnakes1()
                    << " ground truth snakes loaded." << std::endl;
          int desired_rnear = 0;
          int desired_rfar = 0;
          const int max_r = 20;
          double desired_std = 1e6;
          double desired_snr = 0.0;
          double desired_minimum = 0.0;
          double desired_maximum = 0.0;
          for (int rnear = 3; rnear <= max_r; rnear++) {
            for (int rfar = rnear + 1; rfar <= max_r; rfar++) {
              soax::DataContainer snrs;
              multisnake.ComputeGroundTruthLocalSNRs(rnear, rfar, snrs);
              double snr = soax::Mean(snrs);
              double std = soax::StandardDeviation(snrs, snr);
              if (std < desired_std) {
                desired_std = std;
                desired_snr = snr;
                desired_rnear = rnear;
                desired_rfar = rfar;
                desired_minimum = soax::Minimum(snrs);
                desired_maximum = soax::Maximum(snrs);
              }
            }
          }
          outfile << it->filename() << "," << desired_rnear << "," << desired_rfar
                  << "," << desired_snr << "," << desired_std << ","
                  << desired_minimum << "," << desired_maximum << std::endl;
        }
      }
      outfile.close();
    } else {
      std::cout << image_dir << " does not exist." << std::endl;
    }
  } catch (const fs::filesystem_error &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
