/*
 * File: batch_soax.cc
 *
 * This file implements the batch processing of SOAX commandline
 * program.  The input images can be more than one, and the parameters
 * can vary.
 *
 * Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University.
 */

#include <sstream>
#include <iomanip>
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
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

    po::options_description required("Required options");
    required.add_options()
        ("image,i", po::value<std::string>()->required(),
         "Directory or path of input image files")
        ("parameter,p",
         po::value<std::string>()->required(),
         "Path of default parameter file")
        ("snake,s", po::value<std::string>()->required(),
         "Directory or path of output snake files");

    soax::DataContainer ridge_range, stretch_range;
    po::options_description optional("Optional options");
    optional.add_options()
        ("ridge",
         po::value<soax::DataContainer>(&ridge_range)->multitoken(),
         "Range of ridge threshold (start step end)")
        ("stretch",
         po::value<soax::DataContainer>(&stretch_range)->multitoken(),
         "Range of stretching factor (start step end)")
        ("invert", "Use inverted image intensity");

    po::options_description all("Allowed options");
    all.add(generic).add(required).add(optional);
    po::variables_map vm;
    po::store(parse_command_line(argc, argv, all), vm);

    if (vm.count("version")) {
      const std::string version_msg(
          "SOAX Batch 3.5.3\n"
          "Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University.");
      std::cout << version_msg << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("help")) {
      std::cout << "SOAX Commandline Mode: \n" << all;
      return EXIT_SUCCESS;
    }

    po::notify(vm);

    namespace fs = boost::filesystem;
    fs::path image_path(vm["image"].as<std::string>());
    if (!fs::exists(image_path)) {
      std::cerr << image_path << " does not exist. Abort." << std::endl;
      return EXIT_FAILURE;
    }

    fs::path parameter_path(vm["parameter"].as<std::string>());
    if (!fs::exists(parameter_path)) {
      std::cerr << parameter_path << " does not exist. Abort." << std::endl;
      return EXIT_FAILURE;
    }

    fs::path snake_path(vm["snake"].as<std::string>());
    if (!fs::exists(snake_path)) {
      std::cerr << snake_path << " does not exist. Abort." << std::endl;
      return EXIT_FAILURE;
    }

    try {
      soax::Multisnake multisnake;
      if (vm.count("ridge") && vm.count("stretch")) {
        std::cout << "Varying ridge threshold and stretch factor." << std::endl;

        if (fs::is_regular_file(image_path)) {
          std::cout << "Input is single image" << std::endl;

          multisnake.LoadImage(image_path.string());
          if (vm.count("invert"))  multisnake.InvertImageIntensity();
          multisnake.LoadParameters(parameter_path.string());
          multisnake.ComputeImageGradient();

          double ridge_threshold = ridge_range[0];
          while (ridge_threshold < ridge_range[2]) {
            // for (int ridge_exp = ridge_range[0]; ridge_exp < ridge_range[2];
            //      ridge_exp += static_cast<int>(ridge_range[1])) {
            //   double ridge_threshold = std::pow(2, ridge_exp);
            std::cout << "\nSegmentation started on " << image_path
                      << "\nridge_threshold is set to: " << ridge_threshold
                      << std::endl;
            multisnake.set_ridge_threshold(ridge_threshold);
            double stretch = stretch_range[0];
            while (stretch < stretch_range[2]) {
              // for (int stretch_exp = stretch_range[0];
              //      stretch_exp < stretch_range[2];
              //      stretch_exp += static_cast<int>(stretch_range[1])) {
              //   double stretch = std::pow(2, stretch_exp);
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
                  image_path.string(), ridge_threshold, stretch);
              multisnake.SaveSnakes(multisnake.converged_snakes(),
                                    snake_path.string() + snake_name);

              std::cout << "Segmentation completed (Evolution time: "
                        << time_elasped << "s)" << std::endl;
              multisnake.junctions().Reset();
              stretch += stretch_range[1];
            }
            ridge_threshold += ridge_range[1];
          }

        } else if (fs::is_directory(image_path)) {
          std::cout << "Input may contain multiple images." << std::endl;

          typedef std::vector<fs::path> Paths;
          Paths image_paths;
          std::copy(fs::directory_iterator(image_path), fs::directory_iterator(),
                    back_inserter(image_paths));
          std::sort(image_paths.begin(), image_paths.end());
          for (Paths::const_iterator image_it(image_paths.begin());
               image_it != image_paths.end(); ++image_it) {
            // fs::directory_iterator image_end_it;
            // for (fs::directory_iterator image_it(image_path);
            //      image_it != image_end_it; ++image_it) {
            std::string suffix = GetImageSuffix(image_it->string());
            if (suffix != "mha" && suffix != "tif" && suffix != "TIF") {
              std::cout << "Unknown image type: " << suffix << std::endl;
              continue;
            }
            multisnake.LoadImage(image_it->string());
            if (vm.count("invert"))  multisnake.InvertImageIntensity();
            multisnake.LoadParameters(parameter_path.string());
            multisnake.ComputeImageGradient();
            // vary ridge_threshold and stretch
            double ridge_threshold = ridge_range[0];
            while (ridge_threshold < ridge_range[2]) {
              std::cout << "\nSegmentation started on " << *image_it
                        << "\nridge_threshold is set to: " << ridge_threshold
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
                    image_it->string(), ridge_threshold, stretch);
                multisnake.SaveSnakes(multisnake.converged_snakes(),
                                      snake_path.string() + snake_name);

                std::cout << "Segmentation completed (Evolution time: "
                          << time_elasped << "s)" << std::endl;
                multisnake.junctions().Reset();
                stretch += stretch_range[1];
              }
              ridge_threshold += ridge_range[1];
            }
          }
        } else {
          std::cout << image_path
                    << " is neither a regular file nor a directory"
                    << std::endl;
        }
      } else if (!vm.count("ridge") && !vm.count("stretch")) {
        std::cout << "Fixed parameters." << std::endl;
        if (fs::is_regular_file(image_path)) {
          std::cout << "Input is a single image. Please use the GUI version."
                    << std::endl;
        } else if (fs::is_directory(image_path)) {
          std::cout << "Input may contain multiple images." << std::endl;

          typedef std::vector<fs::path> Paths;
          Paths image_paths;
          std::copy(fs::directory_iterator(image_path),
                    fs::directory_iterator(),
                    back_inserter(image_paths));
          std::sort(image_paths.begin(), image_paths.end());
          for (Paths::const_iterator image_it(image_paths.begin());
               image_it != image_paths.end(); ++image_it) {
            // fs::directory_iterator image_end_it;
            // for (fs::directory_iterator image_it(image_path);
            //      image_it != image_end_it; ++image_it) {
            std::string suffix = GetImageSuffix(image_it->string());
            if (suffix != "mha" && suffix != "tif" && suffix != "TIF") {
              std::cout << "Unknown image type: " << suffix << std::endl;
              continue;
            }
            multisnake.LoadImage(image_it->string());
            if (vm.count("invert"))  multisnake.InvertImageIntensity();
            multisnake.LoadParameters(parameter_path.string());
            multisnake.ComputeImageGradient();

            std::cout << "\nSegmentation started on " << *image_it << std::endl;
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

            std::string path_str = image_it->string();
            std::string::size_type slash_pos = path_str.find_last_of("/\\");
            std::string::size_type dot_pos = path_str.find_last_of(".");
            std::string extracted_name = path_str.substr(
                slash_pos+1, dot_pos-slash_pos-1);
            multisnake.SaveSnakes(multisnake.converged_snakes(),
                                  snake_path.string() + extracted_name + ".txt");

            std::cout << "Segmentation completed (Evolution time: "
                      << time_elasped << "s)" << std::endl;
            multisnake.junctions().Reset();
          }
        } else {
          std::cout << image_path << " exists, but is neither a regular file nor a directory" << std::endl;
        }
      } else {
        std::cerr << "Cannot varing one parameter only." << std::endl;
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
  std::ostringstream buffer;
  // std::cout << "default: " << buffer.precision() << std::endl;
  buffer.precision(2);
  // std::cout << "current: " << buffer.precision() << std::endl;
  buffer << std::showpoint << extracted_name << "--ridge"
         << ridge_threshold << "--stretch" << stretch << ".txt";
  return buffer.str();
}


std::string GetImageSuffix(const std::string &image_path) {
  std::string::size_type dot_pos = image_path.find_last_of(".");
  // std::cout << image_path.substr(dot_pos+1) << std::endl;
  return image_path.substr(dot_pos+1);
}
