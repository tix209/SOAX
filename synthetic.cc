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
    std::cout << "Usage: ./syn <original-image> <snake> <output-dir>"
              << std::endl;
    return -1;
  }

  soax::Multisnake multisnake;
  multisnake.LoadImage(argv[1]);
  multisnake.LoadGroundTruthSnakes(argv[2]);
  std::cout << multisnake.GetNumberOfComparingSnakes1()
            << " ground truth snakes loaded." << std::endl;

  const unsigned background = 200;

  // double sigma = 10.0;
  // unsigned foreground = 40;
  // std::ostringstream buffer;
  // buffer << "fg" << foreground << "-sigma" << sigma << ".mha";
  // multisnake.GenerateSyntheticImage(foreground, background, sigma,
  //                                   std::string(argv[3]) + buffer.str());

  for (int i = 2; i < 11; i+=2) {
    double sigma = static_cast<double>(i);

    for (unsigned i = 1; i < 6; i++) {
      unsigned foreground = 20 * i;
      std::ostringstream buffer;
      buffer << "fg" << foreground << "-sigma" << sigma << ".mha";
      multisnake.GenerateSyntheticImage(foreground, background, sigma,
          std::string(argv[3]) + buffer.str());
    }
  }

  // const unsigned foreground = 100;
  // const double sigma = 5.5;
  // std::ostringstream buffer;
  // buffer << "fg" << foreground << "-sigma" << sigma << ".mha";

  // multisnake.GenerateSyntheticImage(foreground, background, sigma,
  //     std::string(argv[3]) + buffer.str());

  return 0;
}
