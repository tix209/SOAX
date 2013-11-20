#include <iostream>
#include <iomanip>
#include <fstream>
#include "boost/program_options.hpp"
#include "multisnake.h"
#include "utility.h"




/*
 * Evaluate the resultant snakes using Vertex Error and Hausdorff
 * Distance if the ground truths are known and using F-function
 * otherwise.
 */

int main(int argc, char **argv) {
  try {
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
        ("version,v", "Print version and exit")
        ("help,h", "Print help message and exit")
        ;
    std::string snake_path, function_path;
    po::options_description required("Required options");
    required.add_options()
        ("snake,s", po::value<std::string>(&snake_path)->required(),
         "Snake file path")
        ("function,f",
         po::value<std::string>(&function_path)->required(),
         "Evaluation function output file path")
        ;

    std::string comparing_snake_path;
    double snr_threshold(0.0);
    double penalizer(0.0);

    po::options_description optional("Optional options");
    optional.add_options()
        ("comparing,c",
         po::value<std::string>(&comparing_snake_path)->required(),
         "Comparing snake file path")
        ("snr,t", po::value<double>(&snr_threshold)->required(),
         "Low SNR threshold")
        ("penalizer,p", po::value<double>(&penalizer)->required(),
         "Constant c for penalizing the low SNR snake points.")
        // ("grad-diff", "Set grad-diff as a variable")
        // ("stretch", "Set stretch as a variable")
        ;

    po::options_description all("Allowed options");
    all.add(generic).add(required).add(optional);

    po::variables_map vm;
    po::store(parse_command_line(argc, argv, all), vm);

    if (vm.count("version")) {
      std::cout << "SOAX Evaluation 1.3\n"
          "Copyright (C) 2013 Ting Xu, IDEA Lab, Lehigh University."
                << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("help")) {
      std::cout << "Usage for SOAX Evaluation: \n";
      std::cout << all;
      return EXIT_SUCCESS;
    }

    po::notify(vm);
    soax::Multisnake multisnake;
    multisnake.LoadImage(soax::GetImagePath(snake_path));
    multisnake.LoadConvergedSnakes(snake_path);

    if (vm.count("comparing")) {
      multisnake.LoadGroundTruthSnakes(comparing_snake_path);
      multisnake.EvaluateByVertexErrorHausdorffDistance(snake_path,
                                                        function_path);
    } else {
      multisnake.EvaluateByFFunction(snr_threshold, penalizer,
                                     snake_path, function_path);
    }
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
