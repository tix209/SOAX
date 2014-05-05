/**
 * File: best_snake.cc
 *
 * This program finds the best resultant snake (evaluated by the
 * F-funtion) under different values of low SNR threshold t and
 * penalizing factor c. If a ground truth snake is present, it also
 * computes the Vertex Error and Hausdorff distance of these best
 * snakes.
 *
 * Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University.
 */


#include <iostream>
#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <fstream>
// #include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "multisnake.h"
#include "utility.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

typedef std::vector<fs::path> Paths;
typedef std::map<std::pair<double, double>, double> TCFMap;
typedef std::map<std::string, TCFMap> FValuesMap;
typedef std::map<std::pair<double, double>, std::string> TCFilenameMap;
typedef std::map<std::pair<double, double>,
                 std::pair<double, double> > TCErrorMap;
typedef std::set<std::string> StringSet;


void CreateSortedSnakePaths(const fs::path &snake_dir, Paths &snake_paths);
void ComputeFilenameFValuesMap(soax::Multisnake &ms,
                               const Paths &snake_paths,
                               int rnear, int rfar,
                               const std::vector<double> &t_range,
                               const std::vector<double> &c_range,
                               FValuesMap &fmap);
void ComputeBestSnakes(soax::Multisnake &ms,
                       const FValuesMap &fmap,
                       int rnear, int rfar,
                       const std::vector<double> &t_range,
                       const std::vector<double> &c_range,
                       bool gt,
                       TCFilenameMap &tc_filenames,
                       TCErrorMap &tc_errors,
                       StringSet &candidate_filenames);
void PrintFMap(const FValuesMap &fmap);
void PrintTCFMap(const TCFMap &m);
void PrintTCFilenameMap(const TCFilenameMap &m, std::ostream &os);
std::string GetBestFilename(const FValuesMap &m, double t, double c,
                            double &min_f);
void PrintTCErrorMap(const TCErrorMap &m, std::ostream &os);


int main (int argc, char **argv) {
  try {
    po::options_description generic("Generic options");
    generic.add_options()
        ("version,v", "Print version and exit")
        ("help,h", "Print help and exit");

    po::options_description required("Required options");
    required.add_options()
        ("image,i", po::value<std::string>()->required(),
         "Path of input image")
        ("snake,s", po::value<std::string>()->required(),
         "Directory of resultant snake files")
        ("output,o", po::value<std::string>()->required(),
         "Path of output 'tc-snake' file")
        ("candidates,a", po::value<std::string>()->required(),
         "Path of output candidate filenames file");

    std::vector<double> t_range(3, 0.0), c_range(3, 0.0);
    t_range[0] = 1.0; // start
    t_range[1] = 0.1; // step
    t_range[2] = 3.01; // end
    c_range[0] = 1.0; // start
    c_range[1] = 0.1; // step
    c_range[2] = 3.01; // end

    po::options_description optional("Optional options");
    optional.add_options()
        ("rnear,n", po::value<int>()->default_value(3),
         "Inner radius of local background annulus")
        ("rfar,f", po::value<int>()->default_value(6),
         "Outer radius of local background annulus")
        ("t-range,t",
         po::value<std::vector<double> >(&t_range)->multitoken(),
         "Range of low SNR threshold (start step end)"
         " Default: 1.0 0.1 3.01")
        ("c-range,c",
         po::value<std::vector<double> >(&c_range)->multitoken(),
         "Range of penalizing factor (start step end)"
         " Default: 1.0 0.1 3.01")
        ("ground-truth,g", po::value<std::string>(),
         "Path of the ground truth snake")
        ("error,e", po::value<std::string>(),
         "Path of the output 'tc-error' file");

    po::options_description all("Allowed options");
    all.add(generic).add(required).add(optional);
    po::variables_map vm;
    po::store(parse_command_line(argc, argv, all), vm);

    if (vm.count("version")) {
      const std::string version_msg(
          "best_snake 1.0\n"
          "Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University.");
      std::cout << version_msg << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("help")) {
      std::cout << "Usage for best_snake: \n" << all;
      return EXIT_SUCCESS;
    }

    po::notify(vm);

    if (vm.count("ground-truth") && !vm.count("error")) {
      std::cerr << "Path of output 'tc-error' file is needed"
          " when ground truth is provided." << std::endl;
      return EXIT_FAILURE;
    }

    try {
      fs::path snake_dir(vm["snake"].as<std::string>());
      if (fs::exists(snake_dir)) {
        soax::Multisnake ms;
        ms.LoadImage(vm["image"].as<std::string>());
        Paths snake_paths;
        CreateSortedSnakePaths(snake_dir, snake_paths);
        FValuesMap fmap;
        int rnear = vm["rnear"].as<int>();
        int rfar = vm["rfar"].as<int>();
        std::cout << rnear << std::endl;
        std::cout << rfar << std::endl;
        std::cout << "t-range: " << t_range[0] << ' ' << t_range[1]
                  << ' ' << t_range[2] << std::endl;
        std::cout << "c-range: " << c_range[0] << ' ' << c_range[1]
                  << ' ' << c_range[2] << std::endl;
        ComputeFilenameFValuesMap(ms, snake_paths, rnear, rfar,
                                  t_range, c_range, fmap);

        TCFilenameMap tc_filenames;
        TCErrorMap tc_errors;
        // std::string gt_path;
        if (vm.count("ground-truth")) {
          ms.LoadGroundTruthSnakes(vm["ground-truth"].as<std::string>());
          std::cout << ms.GetNumberOfComparingSnakes1()
                    << " ground truth snakes loaded." << std::endl;
        }

        StringSet candidate_filenames;

        ComputeBestSnakes(ms, fmap, rnear, rfar, t_range, c_range,
                          vm.count("ground-truth"), tc_filenames,
                          tc_errors, candidate_filenames);

        // PrintFMap(fmap);
        std::ofstream tc_filenames_file(
            vm["output"].as<std::string>().c_str());
        PrintTCFilenameMap(tc_filenames, tc_filenames_file);

        if (vm.count("error")) {
          std::ofstream tc_error_file(vm["error"].as<std::string>().c_str());
          PrintTCErrorMap(tc_errors, tc_error_file);
        }

        std::ofstream candidate_file(
            vm["candidates"].as<std::string>().c_str());
        for (StringSet::const_iterator it = candidate_filenames.begin();
             it != candidate_filenames.end(); ++it) {
          candidate_file << *it << std::endl;
        }
      } else {
        std::cout << snake_dir << " does not exist." << std::endl;
      }
    } catch (const fs::filesystem_error &e) {
      std::cout << e.what() << std::endl;
      return EXIT_FAILURE;
    }
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}

void CreateSortedSnakePaths(const fs::path &snake_dir,
                            Paths &snake_paths) {
  std::copy(fs::directory_iterator(snake_dir),
            fs::directory_iterator(),
            back_inserter(snake_paths));
  std::sort(snake_paths.begin(), snake_paths.end());
}

void ComputeFilenameFValuesMap(soax::Multisnake &ms,
                               const Paths &snake_paths,
                               int rnear, int rfar,
                               const std::vector<double> &t_range,
                               const std::vector<double> &c_range,
                               FValuesMap &fmap) {
  for (Paths::const_iterator it = snake_paths.begin();
       it != snake_paths.end(); ++it) {
    // std::cout << it->filename() << std::endl;
    ms.LoadConvergedSnakes(it->string());
    // ms.PrintSnakes(ms.converged_snakes());
    soax::DataContainer snrs;
    ms.ComputeResultSnakesLocalSNRs(rnear, rfar, snrs);
    // std::cout << rnear << ", " << rfar << ", "
    //           << snrs[0] << std::endl;
    // std::cout << ms.GetNumberOfConvergedSnakes() << std::endl;
    TCFMap tcf;
    for (double t = t_range[0]; t < t_range[2]; t += t_range[1]) {
      for (double c = c_range[0]; c < c_range[2]; c += c_range[1]) {
        double fvalue = ms.ComputeFValue(snrs, t, c);
        // std::cout << "[" << t << ", " << c << "]: " << fvalue
        //           << "\tsnrs: " << soax::Mean(snrs) << std::endl;
        tcf[std::make_pair(t, c)] = fvalue;
      }
    }
    fmap[it->string()] = tcf;
  }
}

void ComputeBestSnakes(soax::Multisnake &ms,
                       const FValuesMap &fmap,
                       int rnear, int rfar,
                       const std::vector<double> &t_range,
                       const std::vector<double> &c_range,
                       bool gt,
                       TCFilenameMap &tc_filenames,
                       TCErrorMap &tc_errors,
                       StringSet &candidate_filenames) {
  for (double t = t_range[0]; t < t_range[2]; t += t_range[1]) {
    for (double c = c_range[0]; c < c_range[2]; c += c_range[1]) {
      double min_f = 1e8;
      std::string filename = GetBestFilename(fmap, t, c, min_f);
      tc_filenames[std::make_pair(t, c)] = filename;

      if ((t + c > 3.0) && (t + c < 5.0))
        candidate_filenames.insert(filename);

      if (gt) {
        ms.LoadConvergedSnakes(filename);
        double vertex_error(100.0), hausdorff(100.0);
        ms.ComputeResultSnakesVertexErrorHausdorffDistance(vertex_error,
                                                           hausdorff);
        tc_errors[std::make_pair(t, c)] =
            std::make_pair(vertex_error, hausdorff);
      }
    }
  }
}


void PrintFMap(const FValuesMap &fmap) {
  for (FValuesMap::const_iterator it = fmap.begin();
       it != fmap.end(); ++it) {
    std::cout << it->first << ":\n";
    // PrintVector(it->second);
    PrintTCFMap(it->second);
    std::cout << std::endl;
  }
}

void PrintTCFMap(const TCFMap &m) {
  for (TCFMap::const_iterator it = m.begin();
       it != m.end(); ++it) {
    std::cout << "(" << it->first.first << ", " << it->first.second
              << "): " << it->second << '\t';
  }
}
// void PrintVector(const std::vector<double> &v) {
//   std::cout << "[ ";
//   for (int i = 0; i < v.size(); i++) {
//     std::cout << v[i] << " ";
//   }
//   std::cout << "]" << std::endl;
// }

// std::string GetMinimalErrorFilename(const FValuesMap &m, unsigned index) {
//   double min_f = 1e8;
//   std::string filename;
//   for (FValuesMap::const_iterator it = m.begin();
//        it != m.end(); ++it) {
//     if (it->second[index] < min_f) {
//       min_f = it->second[index];
//       filename = it->first;
//     }
//   }
//   // std::cout << filename << std::endl;
//   return filename;
// }
void PrintTCFilenameMap(const TCFilenameMap &m, std::ostream &os) {
  for (TCFilenameMap::const_iterator it = m.begin();
       it != m.end(); ++it) {
    os << it->first.first << "\t" << it->first.second << "\t"
       << it->second << std::endl;
  }
}

std::string GetBestFilename(const FValuesMap &m, double t, double c,
                            double &min_f) {
  std::string filename;
  min_f = 1e8;
  for (FValuesMap::const_iterator it = m.begin();
       it != m.end(); ++it) {
    TCFMap::const_iterator f_it = it->second.find(std::make_pair(t, c));
    if (f_it != it->second.end()) {
      if (f_it->second < min_f) {
        min_f = f_it->second;
        filename = it->first;
      }
    } else {
      std::cerr << "cannot find F-value for (" << t << ", " << c << ")."
                << std::endl;
    }
  }

  return filename;
}

void PrintTCErrorMap(const TCErrorMap &m, std::ostream &os) {
  for (TCErrorMap::const_iterator it = m.begin();
       it != m.end(); ++it) {
    os << it->first.first << "\t" << it->first.second << "\t"
       << it->second.first << "\t" << it->second.second << std::endl;
  }
}
