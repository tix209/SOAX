/*
 * File: batch_soax.cc
 *
 * This file implements the batch processing of SOAX commandline program.
 * The input images can be more than one, and the parameters can vary.
 *
 * Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University.
 */

#include <sstream>
#include "boost/program_options.hpp"
#include <boost/filesystem.hpp>
#include "multisnake.h"

std::string ConstructSnakeFilename(const std::string &image_path,
                                   double ridge_threshold, double stretch);
std::string GetImageSuffix(const std::string &image_path);

int main(int argc, char **argv) {
  try {
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
        ("version,v", "Print version and exit")
        ("help,h", "Print help and exit");

    std::string parameter_path, snakes_dir;
    soax::DataContainer ridge_range, stretch_range;
    po::options_description required("Required options");
    required.add_options()
        ("image,i", po::value<std::string>()->required(),
         "Directory of input image files")
        ("parameter,p",
         po::value<std::string>(&parameter_path)->required(),
         "Path of default parameter file")
        ("snake,s", po::value<std::string>(&snakes_dir)->required(),
         "Directory of output snake files")
        ("ridge",
         po::value<soax::DataContainer>(&ridge_range)->multitoken()->
         required(),
         "Range of ridge threshold for SOAC initialization (start step end)")
        ("stretch",
         po::value<soax::DataContainer>(&stretch_range)->multitoken()->
         required(),
         "Range of stretching factor for SOAC evolution (start step end)");

    // soax::DataContainer ridge_range, stretch_range;
    // po::options_description optional("Optional options");
    // optional.add_options()
    //     ("ridge",
    //      po::value<soax::DataContainer>(&ridge_range)->multitoken(),
    //      "Range of ridge threshold for SOAC initialization (start step end)")
    //     ("stretch",
    //      po::value<soax::DataContainer>(&stretch_range)->multitoken(),
    //      "Range of stretching factor for SOAC evolution (start step end)");

    po::options_description all("Allowed options");
    all.add(generic).add(required);
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

      soax::Multisnake multisnake;
      // vary the image
      fs::directory_iterator image_end_it;
      for (fs::directory_iterator image_it(image_dir);
           image_it != image_end_it; ++image_it) {
        std::string suffix = GetImageSuffix(image_it->path().string());
        if (suffix != "mha") {
          std::cout << "Unknown image type: " << suffix << std::endl;
          continue;
        }
        multisnake.LoadImage(image_it->path().string());
        multisnake.LoadParameters(parameter_path);
        multisnake.ComputeImageGradient();
        // vary ridge_threshold and stretch
        double ridge_threshold = ridge_range[0];
        while (ridge_threshold < ridge_range[2]) {
          std::cout << "\nSegmentation started ..." << std::endl;
          std::cout << "ridge_threshold is set to: " << ridge_threshold
                    << std::endl;
          multisnake.set_ridge_threshold(ridge_threshold);
          double stretch = stretch_range[0];
          while (stretch < stretch_range[2]) {
            std::cout << "stretch is set to: " << stretch << std::endl;
            soax::Snake::set_stretch_factor(stretch);
            std::cout << "=========== Current Parameters ==========="
                      << std::endl;
            multisnake.WriteParameters(std::cout);
            std::cout << "=========================================="
                      << std::endl;
            multisnake.InitializeSnakes();

            time_t start, end;
            time(&start);
            multisnake.DeformSnakes();
            time(&end);
            double time_elasped = difftime(end, start);

            multisnake.CutSnakesAtTJunctions();
            multisnake.GroupSnakes();
            std::string snake_name = ConstructSnakeFilename(
                image_it->path().string(), ridge_threshold, stretch);
            multisnake.SaveSnakes(multisnake.converged_snakes(),
                                  snakes_dir + snake_name);
            std::cout << "Segmentation completed (Evolution time: "
                      << time_elasped << "s)" << std::endl;
            multisnake.junctions().Reset();
            stretch += stretch_range[1];
          }
          ridge_threshold += ridge_range[1];
        }
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


std::string ConstructSnakeFilename(const std::string &image_path,
                                   double ridge_threshold, double stretch) {
  std::string::size_type slash_pos = image_path.find_last_of("/\\");
  std::string::size_type dot_pos = image_path.find_last_of(".");
  std::string extracted_name = image_path.substr(
      slash_pos+1, dot_pos-slash_pos-1);
  // std::cout << "extracted name: " << extracted_name << std::endl;
  std::stringstream buffer;
  buffer << extracted_name << "--ridge" << ridge_threshold
         << "--stretch" << stretch << ".txt";
  return buffer.str();
}


std::string GetImageSuffix(const std::string &image_path) {
  std::string::size_type dot_pos = image_path.find_last_of(".");
  // std::cout << image_path.substr(dot_pos+1) << std::endl;
  return image_path.substr(dot_pos+1);
}
