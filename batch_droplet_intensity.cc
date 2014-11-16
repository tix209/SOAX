/** This file computing the mean droplet intensity for all the
 * droplets in a directory. The meta information of droplets such as
 * center location, radius, etc. are provided by another file. The
 * droplet intensity is substracted by an estimated common mean
 * background intensity value.
 *
 * Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University.
 */


#include <iostream>
#include <fstream>
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "multisnake.h"
#include "droplet_info_dict.h"

int main (int argc, char **argv) {
  try {
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
        ("version,v", "Print version and exit")
        ("help,h", "Print help and exit");

    po::options_description required("Required options");
    required.add_options()
        ("image,i", po::value<std::string>()->required(),
         "Directory of images files")
        ("droplet,d", po::value<std::string>()->required(),
         "Path of the droplet info file (.csv)")
        ("background,b", po::value<int>()->required(),
         "Mean intensity of background");

    po::options_description optional("Optional opitions");
    optional.add_options()
        ("output,o", po::value<std::string>(), "Path of output file (.csv)");

    po::options_description all("Allowed options");
    all.add(generic).add(required).add(optional);
    po::variables_map vm;
    po::store(parse_command_line(argc, argv, all), vm);

    if (vm.count("version")) {
      std::string version_msg("Batch droplet intensity 0.1\n"
                              "Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University.");
      std::cout << version_msg << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("help")) {
      std::cout << all;
      return EXIT_SUCCESS;
    }

    po::notify(vm);

    namespace fs = boost::filesystem;
    fs::path image_dir(vm["image"].as<std::string>());
    if (!fs::exists(image_dir)) {
      std::cerr << image_dir << " does not exist. Abort." << std::endl;
      return EXIT_FAILURE;
    }

    try {
      soax::Multisnake multisnake;
      if (fs::is_directory(image_dir)) {
        soax::DropletInfoDict droplet_infos;
        std::string info_filename(vm["droplet"].as<std::string>());
        if (!droplet_infos.LoadInfo(info_filename)) {
          std::cerr << "Error reading " << info_filename << std::endl;
          return EXIT_FAILURE;
        }

        std::string output_filename("droplet_mean_foregrounds.txt");
        if (vm.count("output"))
          output_filename = vm["output"].as<std::string>();
        fs::path outpath(output_filename);
        if (fs::exists(outpath)) {
          std::cout << "Warning: output file aleady exists. "
              "New output will be appended to its end!" << std::endl;
        }

        std::ofstream outfile(output_filename.c_str(), std::ofstream::app);
        if (!outfile) {
          std::cerr << "Couldn't open outfile " << output_filename << std::endl;
          return EXIT_FAILURE;
        }
        outfile << "ImageFileName,MeanForeground" << std::endl;

        typedef std::vector<fs::path> Paths;
        Paths image_paths;
        std::copy(fs::directory_iterator(image_dir), fs::directory_iterator(),
                  back_inserter(image_paths));
        std::sort(image_paths.begin(), image_paths.end());
        for (Paths::const_iterator it = image_paths.begin();
             it != image_paths.end(); ++it) {
          std::string imagename = it->filename().string();
          if (droplet_infos.Has(imagename)) {
            multisnake.LoadImage(it->string());
            soax::PointType center = droplet_infos.Get(imagename).center();
            double radius = droplet_infos.Get(imagename).radius();
            double intensity = multisnake.ComputeDropletMeanIntensity(center, radius);
            int background = vm["background"].as<int>();
            outfile << imagename << "," << intensity - background << std::endl;
          }
        }
        outfile.close();
      } else {
        std::cerr << "Argument of --image is not a directory. Abort." << std::endl;
        return EXIT_FAILURE;
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
