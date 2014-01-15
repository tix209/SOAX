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
              << std::endl;
    return EXIT_FAILURE;
  }

  namespace fs = boost::filesystem;
  std::string image_path = argv[1];
  std::string gt_snake_path = argv[2];
  // std::string snake_dir = argv[3];
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

      for (int i = 1; i < 10; ++i) {
        double threshold = 0.1 * i * snr;
        for (int j = 0; j < 10; ++j) {
          double penalizer = 1.0 + 0.5 * j;
          double gt_fvalue = multisnake.ComputeGroundTruthFValue(threshold,
                                                                 penalizer,
                                                                 4, 12);
          bool snakes_greater_fvalue = true;

          fs::directory_iterator end_it;
          for (fs::directory_iterator it(snake_dir); it != end_it; ++it) {
            // std::cout << "resultant snake: " << it->path().filename()
            //           << std::endl;
            multisnake.LoadConvergedSnakes(it->path().string());
            // std::cout << multisnake.GetNumberOfConvergedSnakes()
            //           << " resultant snakes loaded." << std::endl;

            double fvalue = multisnake.ComputeResultSnakesFValue(threshold,
                                                                 penalizer,
                                                                 4, 12);
            if (fvalue < gt_fvalue) {
              snakes_greater_fvalue = false;
              // std::cout << "t: " << i*0.1 << "\t" << "c: "
              //           << penalizer << "\tgt fvalue: " << gt_fvalue
              //           << "\tresult fvalue: " << fvalue << std::endl;
              // std::cout << it->path().filename() << std::endl;
              break;
            }
          }

          if (snakes_greater_fvalue) {
            std::cout << "t: " << i*0.1 << "\t"
                      << "c: " << penalizer << std::endl;
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
