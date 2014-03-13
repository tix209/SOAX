/*
 * This file implements the batch generation of snapshots of a vtk
 * rendering window without showing it.
 */

#include <iostream>
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "multisnake.h"
#include "viewer.h"


int main(int argc, char **argv) {
  try {
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
        ("version,v", "Print version and exit")
        ("help,h", "Print help and exit");

    std::string snake_dir, output_dir;
    po::options_description required("Required options");
    required.add_options()
        ("image,i", po::value<std::string>()->required(), "Directory of image files")
        ("snake,s", po::value<std::string>(&snake_dir)->required(), "Directory of snake files")
        ("output,o", po::value<std::string>(&output_dir)->required(), "Output directory");
    po::options_description all("Allowed options");
    all.add(generic).add(required);
    po::variables_map vm;
    po::store(parse_command_line(argc, argv, all), vm);
    if (vm.count("version")) {
      const std::string version_msg(
          "Batch snapshots 1.0\n"
          "Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University.");
      std::cout << version_msg << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("help")) {
      std::cout << "Batch snapshots usage: \n" << all;
      return EXIT_SUCCESS;
    }
    po::notify(vm);

    namespace fs = boost::filesystem;
    fs::path image_path(vm["image"].as<std::string>());
    if (!fs::exists(image_path)) {
      std::cerr << image_path << " does not exist. Abort." << std::endl;
      return EXIT_FAILURE;
    }

    // fs::path snake_path(vm["snake"].as<std::string>());
    // if (!fs::exists(snake_path)) {
    //   std::cerr << image_path << " does not exist. Abort." << std::endl;
    //   return EXIT_FAILURE;
    // }

    // fs::path output_path(vm["output"].as<std::string>());
    // if (!fs::exists(output_path)) {
    //   std::cerr << output_path << " does not exist. Abort." << std::endl;
    //   return EXIT_FAILURE;
    // }

    QApplication app(argc, argv);
    soax::MainWindow window;
    soax::Multisnake multisnake;
    soax::Viewer viewer;

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
      multisnake.LoadImage(image_it->string());
      viewer.SetupImage(multisnake.image());
      viewer.ToggleSlicePlanes(true);
      viewer.ToggleMIPRendering(false);
      viewer.ToggleOrientationMarker(false);
      viewer.ToggleCornerText(false);
      viewer.ToggleBoundingBox(false);
      viewer.ToggleCubeAxes(false);

      std::string snake_path = snake_dir + multisnake.GetImageName(false) + ".txt";
      std::cout << snake_path << std::endl;
    }
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}
