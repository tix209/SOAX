/*
 * File: evaluate_f.cc
 *
 * This file implements the evaluation of segmentation using
 * F-function. The proper meta-parameters t and c are determined by
 * another file result_f_value.cc.
 *
 * Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University
 */


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
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

    std::string image_path, snake_dir, output_path;
    double t(0.0);
    double c(0.0);
    po::options_description required("Required options");
    required.add_options()
        ("image,i", po::value<std::string>(&image_path)->required(),
         "Image file path")
        ("snake,s", po::value<std::string>(&snake_dir)->required(),
         "Snake file directory")
        ("threshold,t", po::value<double>(&t)->required(),
         "Constant factor for low SNR threshold")
        ("penalizer,c", po::value<double>(&c)->required(),
         "Constant factor for penalizing low SNR snake points.")
        ("output,o", po::value<std::string>(&output_path)->required(),
         "Evaluation output file path")
        ;

    int radial_near(4);
    int radial_far(12);
    po::options_description optional("Optional options");
    optional.add_options()
        ("near,n", po::value<int>(&radial_near),
         "Inner radius of annulus for local SNR estimation. (default: 4)")
        ("far,f", po::value<int>(&radial_far),
         "Outer radius of annulus for local SNR estimation. (default: 12)")
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


    // application code
    namespace fs = boost::filesystem;
    fs::path snakes_path(snake_dir);

    try {
      if (fs::exists(snakes_path)) {
        soax::Multisnake multisnake;
        multisnake.LoadImage(image_path);
        double snr = multisnake.ComputeImageSNR2();
        std::cout << "Image SNR: " << snr << std::endl;
        double threshold = t * snr;

        std::ofstream outfile;
        outfile.open(output_path.c_str());
        if (!outfile.is_open()) {
          std::cerr << "Couldn't open output f-function file: "
                    << output_path << std::endl;
          return EXIT_FAILURE;
        }
        const unsigned width = 15;
        outfile << "Snakes directory\t" << snake_dir << "\n"
                << "Otsu's image SNR\t" << snr << "\n"
                << "Low SNR factor\t" << t << "\n"
                << "Low SNR threshold\t" << threshold << "\n"
                << "Penalizing factor\t" << c << "\n"
                << "Radial near\t" << radial_near << "\n"
                << "Radial far\t" << radial_far << "\n"
                << std::setw(width) << "Ridge"
                << std::setw(width) << "Stretch"
                << std::setw(width) << "Fvalue" << std::endl;

        typedef std::vector<fs::path> Paths;
        Paths sorted_snakes_path;
        std::copy(fs::directory_iterator(snakes_path), fs::directory_iterator(),
                  back_inserter(sorted_snakes_path));
        std::sort(sorted_snakes_path.begin(), sorted_snakes_path.end());

        for (Paths::const_iterator it = sorted_snakes_path.begin();
             it != sorted_snakes_path.end(); ++it) {
          // std::cout << "snake: " << it->path().filename() << std::endl;
          multisnake.LoadConvergedSnakes(it->string());
          // std::cout << multisnake.GetNumberOfConvergedSnakes()
          //           << " resultant snakes loaded." << std::endl;

          // double fvalue = multisnake.ComputeResultSnakesFValue(threshold,
          //                                                      c,
          //                                                      radial_near,
          //                                                      radial_far);
          soax::DataContainer snrs;
          multisnake.ComputeResultSnakesLocalSNRs(radial_near,
                                                  radial_far,
                                                  snrs);
          double fvalue = multisnake.ComputeFValue(snrs, threshold, c);
          outfile << std::setw(width) << multisnake.ridge_threshold()
                  << std::setw(width) << soax::Snake::stretch_factor()
                  << std::setw(width) << fvalue
                  << std::endl;
        }
        outfile.close();
      } else {
        std::cout << snakes_path << " does not exist." << std::endl;
      }
    } catch (const fs::filesystem_error &e) {
      std::cout << e.what() << std::endl;
      return EXIT_FAILURE;
    }
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
