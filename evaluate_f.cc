/*
 * File: evaluate_f.cc
 *
 * This file implements the evaluation of segmentation using F-function.
 *
 * Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University
 */


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "boost/program_options.hpp"
#include "multisnake.h"
#include "utility.h"


int main (int argc, char **argv) {
  try {
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
        ("version,v", "Print version and exit")
        ("help,h", "Print help message and exit")
        ;
    std::string snake_path, output_path;
    po::options_description required("Required options");
    required.add_options()
        ("snake,s", po::value<std::string>(&snake_path)->required(),
         "Snake file path")
        ("output,o", po::value<std::string>(&output_path)->required(),
         "Evaluation output file path")
        ;

    double t(0.0);
    double c(0.0);
    int radial_near(3);
    int radial_far(6);

    po::options_description optional("Optional options");
    optional.add_options()
        ("t,t", po::value<double>(&t), "Constant factor of low SNR threshold")
        ("c,c", po::value<double>(&c),
         "Constant factor for penalizing low SNR snake points.")
        ("near,n", po::value<int>(&radial_near),
         "Inner radius of annulus for local SNR estimation.")
        ("far,f", po::value<int>(&radial_far),
         "Outer radius of annulus for local SNR estimation.")
        ;

    po::options_description all("Allowed options");
    all.add(generic).add(required).add(optional);

    po::variables_map vm;
    po::store(parse_command_line(argc, argv, all), vm);

    if (vm.count("version")) {
      std::cout << "SOAX Evaluation by F-function 1.0\n"
          "Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University."
                << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("help")) {
      std::cout << "Usage for SOAX Evaluation by F-function:\n";
      std::cout << all;
      return EXIT_SUCCESS;
    }

    po::notify(vm);



    // appliation code here
    soax::Multisnake multisnake;
    multisnake.LoadImage(soax::GetImagePath(snake_path));
    multisnake.LoadConvergedSnakes(snake_path);
    // std::cout << multisnake.GetNumberOfConvergedSnakes()
    //           << " resultant snakes loaded." << std::endl;

    double snr = multisnake.ComputeImageSNR();
    std::cout << "Image SNR: " << snr << std::endl;
    double snr_threshold = t * snr;
    multisnake.EvaluateByFFunction(snr_threshold, c, radial_near,
                                   radial_far, snake_path, output_path);
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
