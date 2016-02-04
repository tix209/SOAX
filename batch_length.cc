/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file implements the commandline utility for batch computation of
 * SOACs' length. The input is a directory containing one or more SOAC file;
 * the output is a text file containing the length of all the SOACs in those
 * files.
 */


#include <vector>
#include <fstream>
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "./multisnake.h"

int main(int argc, char **argv) {
  try {
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
        ("version,v", "Print version and exit")
        ("help,h", "Print help and exit");
    po::options_description required("Required options");
    required.add_options()
        ("input,i", po::value<std::string>()->required(),
         "Directory of SOAC files");

    po::options_description optional("Optional opitions");
    optional.add_options()
        ("output,o", po::value<std::string>(), "Output file path");

    po::options_description all("Allowed options");
    all.add(generic).add(required).add(optional);

    po::variables_map vm;
    po::store(parse_command_line(argc, argv, all), vm);

    if (vm.count("version")) {
      std::string version_msg(
          "Batch Length 3.6.0 \n"
          "Computing SOAC lengths (in pixels) from multiple input SOAC files.\n"
          "Copyright (C) 2016, Lehigh University.");
      std::cout << version_msg << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("help")) {
      std::cout << all;
      return EXIT_SUCCESS;
    }
    po::notify(vm);

    namespace fs = boost::filesystem;
    fs::path snake_path(vm["input"].as<std::string>());
    if (!fs::exists(snake_path)) {
      std::cerr << snake_path << " does not exist. Abort." << std::endl;
      return EXIT_FAILURE;
    }
    try {
      std::string output_filename("batch_length_output.csv");
      if (vm.count("output"))
        output_filename = vm["output"].as<std::string>();
      std::ofstream outfile(output_filename.c_str());
      if (!outfile) {
        std::cerr << "Couldn't open " << output_filename << std::endl;
        return EXIT_FAILURE;
      }
      outfile << "SnakeIndex,Length" << std::endl;

      soax::Multisnake multisnake;
      typedef std::vector<fs::path> Paths;
      Paths snake_paths;
      std::copy(fs::directory_iterator(snake_path),
                fs::directory_iterator(),
                back_inserter(snake_paths));
      std::sort(snake_paths.begin(), snake_paths.end());

      for (Paths::const_iterator it(snake_paths.begin());
           it != snake_paths.end(); it++) {
        std::cout << it->string() << std::endl;
        multisnake.LoadConvergedSnakes(it->string());
        std::cout << multisnake.GetNumberOfConvergedSnakes() << std::endl;
        multisnake.ComputeSnakeLength(1.0, outfile);
        multisnake.Reset();
      }
    } catch (const fs::filesystem_error &e) {
      std::cout << e.what() << std::endl;
      return EXIT_FAILURE;
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
