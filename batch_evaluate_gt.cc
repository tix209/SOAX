/*
 * File: batch_evaluate_gt.cc
 *
 * This file implements the evaluation of segmentation using Vertex
 * Error and Hausdorff Distance against the ground truths on multiple
 * images with fixed parameters.
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

    std::string image_dir, gt_snake_dir, result_snake_dir, output_path;
    po::options_description required("Required options");
    required.add_options()
        ("image,i", po::value<std::string>(&image_dir)->required(),
         "Image file directory")
        ("ground-truth,g", po::value<std::string>(&gt_snake_dir)->required(),
         "Ground truth snake directory")
        ("snake,s", po::value<std::string>(&result_snake_dir)->required(),
         "Resultant snake directory")
        ("output,o", po::value<std::string>(&output_path)->required(),
         "Evaluation output file path");

    po::options_description all("Allowed options");
    all.add(generic).add(required);

    po::variables_map vm;
    po::store(parse_command_line(argc, argv, all), vm);

    if (vm.count("version")) {
      std::cout << "Batch SOAX evaluation against ground truth with fixed parameters 1.0\n"
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
        std::ofstream outfile;
        outfile.open(output_path.c_str());
        if (!outfile.is_open()) {
          std::cerr << "Couldn't open output evaluation file: "
                    << output_path << std::endl;
          return EXIT_FAILURE;
        }
        const unsigned width = 20;
        outfile << "Snakes directory\t" << result_snake_dir << "\n"
                << "Ground truth snake\t" << gt_snake_dir << "\n"
                << std::setw(width) << "Image"
                << std::setw(width) << "Vertex"
                << std::setw(width) << "Hausdorff"
                << std::endl;

        typedef std::vector<fs::path> Paths;
        Paths sorted_snakes_path;
        std::copy(fs::directory_iterator(snakes_path), fs::directory_iterator(),
                  back_inserter(sorted_snakes_path));
        // std::sort(sorted_snakes_path.begin(), sorted_snakes_path.end());
        soax::Multisnake multisnake;

        for (Paths::const_iterator it = sorted_snakes_path.begin();
             it != sorted_snakes_path.end(); ++it) {
          // std::cout << "snake: " << it->stem() << std::endl;
          std::string stem = it->stem().string();
          std::string image_filename = stem + ".mha";
          std::cout << image_filename << std::endl;
          multisnake.LoadImage(image_dir + image_filename);
          std::string gt_snake_name = gt_snake_dir +
              stem.substr(0, stem.find('_')) + ".txt";
          std::cout << gt_snake_name << std::endl;
          multisnake.LoadGroundTruthSnakes(gt_snake_name);
          multisnake.LoadConvergedSnakes(it->string());
          std::cout << multisnake.GetNumberOfConvergedSnakes()
                    << " resultant snakes loaded." << std::endl;
          double vertex_error(100.0), hausdorff(100.0);
          multisnake.ComputeResultSnakesVertexErrorHausdorffDistance(
              vertex_error, hausdorff);
          outfile << std::setw(width) << image_filename
                  << std::setw(width) << vertex_error
                  << std::setw(width) << hausdorff << std::endl;
          multisnake.Reset();
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
