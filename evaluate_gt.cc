/*
 * File: evaluate_gt.cc
 *
 * This file implements the evaluation of segmentation using Vertex
 * Error and Hausdorff Distance against the ground truths.
 *
 * Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University
 */


#include <iostream>
#include <iomanip>
#include <fstream>
#include "boost/program_options.hpp"
#include "multisnake.h"
#include "utility.h"


int main(int argc, char **argv) {
  try {
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
        ("version,v", "Print version and exit")
        ("help,h", "Print help message and exit")
        ;
    std::string snake_path, comparing_snake_path, output_path;
    po::options_description required("Required options");
    required.add_options()
        ("snake,s", po::value<std::string>(&snake_path)->required(),
         "Snake file path")
        ("comparing,c", po::value<std::string>(&comparing_snake_path),
         "Comparing snake file path")
        ("output,o", po::value<std::string>(&output_path)->required(),
         "Evaluation output file path")
        ;

    // std::string comparing_snake_path;
    // double snr_threshold(0.0);
    // double penalizer(0.0);
    // int radial_near(3);
    // int radial_far(6);

    po::options_description all("Allowed options");
    all.add(generic).add(required);

    po::variables_map vm;
    po::store(parse_command_line(argc, argv, all), vm);

    if (vm.count("version")) {
      std::cout << "SOAX Evaluation GT 1.3\n"
          "Copyright (C) 2013 Ting Xu, IDEA Lab, Lehigh University."
                << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("help")) {
      std::cout << "Usage for SOAX Evaluation GT: \n";
      std::cout << all;
      return EXIT_SUCCESS;
    }

    po::notify(vm);
    soax::Multisnake multisnake;
    multisnake.LoadImage(soax::GetImagePath(snake_path));
    multisnake.LoadConvergedSnakes(snake_path);
    multisnake.LoadGroundTruthSnakes(comparing_snake_path);
    multisnake.EvaluateByVertexErrorHausdorffDistance(snake_path,
                                                      output_path);
    // multisnake.PrintGroundTruthLocalSNRValues(radial_near, radial_far);
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
