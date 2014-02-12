/*
 * File: result_f_value.cc
 *
 * This file implements finding the suitable t, c meta parameters for
 * F-function evaluation. Proper t, c values ensure that the ground
 * truth snakes gives the minimum f-value among all resultant snakes.
 * Whether these proper values can be found depend on the value of
 * radial_near and radial_far for local SNR estimation. Currently we
 * set radial_near = 4, radial_far = 12, and we use snake point
 * intensity as the foreground intensity.
 *
 * Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University.
 */


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "multisnake.h"


int main (int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "./rsf <image_path> <gt_snake_path> <resultant_snake_dir>"
        " <radial_near> <radial_far>"
              << std::endl;
    return EXIT_FAILURE;
  }

  namespace fs = boost::filesystem;
  std::string image_path = argv[1];
  std::string gt_snake_path = argv[2];
  fs::path snake_dir(argv[3]);

  try {
    if (fs::exists(snake_dir)) {
      soax::Multisnake multisnake;
      multisnake.LoadImage(image_path);
      double snr = multisnake.ComputeImageSNR();
      std::cout << "Image SNR: " << snr << std::endl;
      multisnake.LoadGroundTruthSnakes(gt_snake_path);
      std::cout << multisnake.GetNumberOfComparingSnakes1()
                << " ground truth snakes loaded." << std::endl;

      // unsigned least_number_of_violation = soax::kBigNumber;
      double t = 0.0;
      double c = 0.0;
      int radial_near = atoi(argv[4]);
      int radial_far = atoi(argv[5]);
      std::cout << "radial_near: " << radial_near << "\t"
                << "radial_far: " << radial_far << std::endl;

      soax::DataContainer gt_snrs;
      multisnake.ComputeGroundTruthLocalSNRs(
          radial_near, radial_far, gt_snrs);
      // std::cout << "gt_snr size: " << gt_snr.size() << std::endl;

      std::vector<std::vector<double> > result_snrs_vector;
      fs::directory_iterator end_it;
      // unsigned num_of_snake_sets = 0;
      for (fs::directory_iterator it(snake_dir); it != end_it; ++it) {
        multisnake.LoadConvergedSnakes(it->path().string());
        soax::DataContainer result_snrs;
        multisnake.ComputeResultSnakesLocalSNRs(
            radial_near, radial_far, result_snrs);
        result_snrs_vector.push_back(result_snrs);
        // num_of_snake_sets++;
      }


      for (int i = 1; i <= 10; ++i) {
        double threshold = 0.1 * i * snr;
        // double threshold = 1.5;
        for (int j = 1; j <= 20; ++j) {
          double penalizer = j * 0.5;
          double gt_fvalue = multisnake.ComputeFValue(gt_snrs,
                                                      threshold,
                                                      penalizer);
          bool snakes_greater_fvalue = true;
          // unsigned number_of_violation = 0;
          // std::cout << "t: " << i*0.1 << "\t" << "c: " << penalizer
          //           << "\tgt fvalue: " << gt_fvalue << std::endl;

          fs::directory_iterator end_it;
          unsigned index = 0;
          for (fs::directory_iterator it(snake_dir); it != end_it; ++it) {
            double result_fvalue = multisnake.ComputeFValue(
                result_snrs_vector[index], threshold, penalizer);

            if (result_fvalue < gt_fvalue + soax::kEpsilon) {
              snakes_greater_fvalue = false;
              break;
            }
            index++;
          }

          if (snakes_greater_fvalue) {
            std::cout << "t: " << i*0.1 << "\t";
            std::cout << "c: " << penalizer << std::endl;
          }
        }
      }
    } else {
      std::cout << snake_dir << " does not exist." << std::endl;
    }
  } catch (const fs::filesystem_error &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
