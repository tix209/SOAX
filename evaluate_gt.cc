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
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
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

    std::string image_path, gt_snake_path, result_snake_dir, output_path;
    po::options_description required("Required options");
    required.add_options()
        ("image,i", po::value<std::string>(&image_path)->required(),
         "Image file path")
        ("ground-truth,g", po::value<std::string>(&gt_snake_path)->required(),
         "Ground truth snake file path")
        ("snake,s", po::value<std::string>(&result_snake_dir)->required(),
         "Resultant snake directory")
        ("output,o", po::value<std::string>(&output_path)->required(),
         "Evaluation output file path")
        ;

    po::options_description all("Allowed options");
    all.add(generic).add(required);

    po::variables_map vm;
    po::store(parse_command_line(argc, argv, all), vm);

    if (vm.count("version")) {
      std::cout << "SOAX Evaluation using Ground Truth 1.4\n"
          "Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University."
                << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("help")) {
      std::cout << "Usage: \n";
      std::cout << all;
      return EXIT_SUCCESS;
    }

    po::notify(vm);

    // application code
    namespace fs = boost::filesystem;
    fs::path snakes_path(result_snake_dir);

    try {
      if (fs::exists(snakes_path)) {
        soax::Multisnake multisnake;
        multisnake.LoadImage(image_path);
        multisnake.LoadGroundTruthSnakes(gt_snake_path);
        std::cout << multisnake.GetNumberOfComparingSnakes1()
                  << " ground truth snakes loaded." << std::endl;

        std::ofstream outfile;
        outfile.open(output_path.c_str());
        if (!outfile.is_open()) {
          std::cerr << "Couldn't open output evaluation file: "
                    << output_path << std::endl;
          return EXIT_FAILURE;
        }
        const unsigned width = 15;
        outfile << "Snakes directory\t" << result_snake_dir << "\n"
                << "Ground truth snake\t" << gt_snake_path << "\n"
                << std::setw(width) << "Ridge"
                << std::setw(width) << "Stretch"
                << std::setw(width) << "Vertex"
                << std::setw(width) << "Hausdorff"
                << std::endl;

        fs::directory_iterator end_it;
        for (fs::directory_iterator it(snakes_path); it != end_it; ++it) {
          // std::cout << "snake: " << it->path().filename() << std::endl;
          multisnake.LoadConvergedSnakes(it->path().string());
          // std::cout << multisnake.GetNumberOfConvergedSnakes()
          //           << " resultant snakes loaded." << std::endl;
          double vertex_error(100.0), hausdorff(100.0);
          multisnake.ComputeResultSnakesVertexErrorHausdorffDistance(
              vertex_error, hausdorff);
          outfile << std::setw(width) << multisnake.ridge_threshold()
                  << std::setw(width) << soax::Snake::stretch_factor()
                  << std::setw(width) << vertex_error
                  << std::setw(width) << hausdorff << std::endl;
        }
        outfile.close();
      } else {
        std::cout << snakes_path << " does not exist." << std::endl;
      }
    } catch (const fs::filesystem_error &e) {
      std::cout << e.what() << std::endl;
      return EXIT_FAILURE;
    }
    // multisnake.PrintGroundTruthLocalSNRValues(radial_near, radial_far);
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
