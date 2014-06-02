/*
 * File: synthetic.cc
 *
 * This file implements the generation of synthetic curvilinear
 * network images with uniform foreground and various level of
 * Gaussian noise. The input is a snake file (in JFilament format) and
 * the output are a bunch of images in the specified folder.
 *
 * Copyright (C) Ting Xu 2014, IDEA Lab, Lehigh University
 */

#include <iostream>
#include <sstream>
#include "multisnake.h"

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cout << "Usage: ./syn <original-image-path> <snake-path> "
        "<output-dir>" << std::endl;
    return -1;
  }

  soax::Multisnake multisnake;
  multisnake.LoadImage(argv[1]);
  multisnake.LoadGroundTruthSnakes(argv[2]);
  std::cout << multisnake.GetNumberOfComparingSnakes1()
            << " ground truth snakes loaded." << std::endl;

  const unsigned background = 200;
  const unsigned foreground = 20;

  for (int i = 0; i < 8; i++) {
    double sigma = static_cast<double>(i);
    std::ostringstream buffer;
    buffer << "fg" << foreground << "-bg" << background
           << "-s" << sigma << ".tif";
    std::string output_file_path = std::string(argv[3]) + buffer.str();
    multisnake.GenerateSyntheticImage(foreground, background, sigma,
                                      output_file_path);
    // multisnake.GenerateSyntheticRealImage(fg_ratio, sigma,
    //                                       output_file_path);
  }
  return 0;
}
