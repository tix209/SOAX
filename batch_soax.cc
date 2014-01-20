/*
 * File: batch_soax.cc
 *
 * This file implements the batch processing of SOAX commandline program.
 * The input images can be more than one, and the parameters can vary.
 *
 * Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University.
 */

#include "boost/program_options.hpp"
#include <boost/filesystem.hpp>
#include "multisnake.h"

int main(int argc, char **argv) {
  try {
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
        ("version,v", "Print version and exit")
        ("help,h", "Print help and exit");

    // std::string image_dir, parameter_path, snakes_dir;
    std::string parameter_path;
    po::options_description required("Required options");
    required.add_options()
        ("image,i", po::value<std::string>()->required(),
         "Directory of input image files")
        ("parameter,p",
         po::value<std::string>(&parameter_path)->required(),
         "Path of default parameter file")
        ("snake,s", po::value<std::string>()->required(),
         "Directory of output snake files");

    soax::DataContainer ridge_range, stretch_range;
    po::options_description optional("Optional options");
    optional.add_options()
        ("ridge",
         po::value<soax::DataContainer>(&ridge_range)->multitoken(),
         "Range of ridge threshold for SOAC initialization (start step end)")
        ("stretch",
         po::value<soax::DataContainer>(&stretch_range)->multitoken(),
         "Range of stretching factor for SOAC evolution (start step end)");

    po::options_description all("Allowed options");
    all.add(generic).add(required).add(optional);
    po::variables_map vm;
    po::store(parse_command_line(argc, argv, all), vm);

    if (vm.count("version")) {
      const std::string version_msg(
          "SOAX Batch 3.5.2\n"
          "Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University.");
      std::cout << version_msg << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("help")) {
      std::cout << "SOAX Batch Mode: \n" << all;
      return EXIT_SUCCESS;
    }

    po::notify(vm);

    namespace fs = boost::filesystem;
    fs::path image_dir(vm["image"].as<std::string>());

    try {
      if (!fs::exists(image_dir)) {
        std::cerr << image_dir << " does not exist. Abort." << std::endl;
        return EXIT_FAILURE;
      }

      for (soax::DataContainer::const_iterator it = ridge_range.begin();
           it != ridge_range.end(); ++it)
        std::cout << *it << std::endl;

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
