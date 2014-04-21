#include <iostream>
#include <cstdlib>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "multisnake.h"
#include "utility.h"

typedef std::map<std::string, std::vector<double> > FValuesMap;
typedef std::map<double, std::pair<double, double> > TErrorMap;
std::vector<double> ComputeFValueVector(soax::Multisnake &ms,
                                        double start, double end, double step);
void PrintMap(const FValuesMap &fmap);
void PrintVector(const std::vector<double> &v);
void PrintThresholdErrorMap(const TErrorMap &m);
std::string GetMinimalErrorFilename(const FValuesMap &m, unsigned index);

int main (int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "./bestt <image_path> <gt_path> <snake_dir> " << std::endl;
    return EXIT_FAILURE;
  }

  namespace fs = boost::filesystem;
  std::string image_path = argv[1];
  std::string gt_path = argv[2];
  fs::path snake_dir(argv[3]);

  try {
    if (fs::exists(snake_dir)) {
      soax::Multisnake ms;
      ms.LoadImage(image_path);

      FValuesMap fmap;
      typedef std::vector<fs::path> Paths;
      Paths sorted_snakes_path;
      std::copy(fs::directory_iterator(snake_dir),
                fs::directory_iterator(),
                back_inserter(sorted_snakes_path));
      std::sort(sorted_snakes_path.begin(), sorted_snakes_path.end());

      double t_start = 1.0;
      double t_end = 9.9;
      double t_step = 1.0;
      for (Paths::const_iterator it = sorted_snakes_path.begin();
           it != sorted_snakes_path.end(); ++it) {
        // std::cout << it->filename() << std::endl;
        ms.LoadConvergedSnakes(it->string());
        fmap[it->string()] = ComputeFValueVector(ms, t_start, t_end, t_step);
      }

      TErrorMap t_error_map; // threshold - Hausdorff distance/Vertex error map
      ms.LoadGroundTruthSnakes(gt_path);
      std::cout << ms.GetNumberOfComparingSnakes1()
                << " ground truth snakes loaded." << std::endl;

      unsigned num_entries = static_cast<unsigned>((t_end - t_start) / t_step)+1;
      for (unsigned i = 0; i < num_entries; i++) {
        std::string filename = GetMinimalErrorFilename(fmap, i);
        double t = t_start + t_step * i;
        std::cout << filename << std::endl;
        ms.LoadConvergedSnakes(filename);
        double vertex_error(100.0), hausdorff(100.0);
        ms.ComputeResultSnakesVertexErrorHausdorffDistance(vertex_error, hausdorff);
        // std::cout << hausdorff << std::endl;
        t_error_map[t] = std::make_pair(vertex_error, hausdorff);
      }

      PrintMap(fmap);
      PrintThresholdErrorMap(t_error_map);
    } else {
      std::cout << snake_dir << " does not exist." << std::endl;
    }
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}


std::vector<double> ComputeFValueVector(soax::Multisnake &ms,
                                        double start, double end, double step) {
  const double c = 2.0;
  std::vector<double> v;
  for (double t = start; t < end; t += step) {
    soax::DataContainer snrs;
    ms.ComputeResultSnakesLocalSNRs(4, 12, snrs);
    v.push_back(ms.ComputeFValue(snrs, t, c));
  }
  return v;
}

void PrintThresholdErrorMap(const TErrorMap &m) {
  for (TErrorMap::const_iterator it = m.begin();
       it != m.end(); ++it) {
    std::cout << it->first << ": [" << it->second.first << ' '
              << it->second.second << "]"<< std::endl;
  }
}

void PrintMap(const FValuesMap &fmap) {
  for (FValuesMap::const_iterator it = fmap.begin();
       it != fmap.end(); ++it) {
    std::cout << it->first << ": ";
    PrintVector(it->second);
  }
}

void PrintVector(const std::vector<double> &v) {
  std::cout << "[ ";
  for (int i = 0; i < v.size(); i++) {
    std::cout << v[i] << " ";
  }
  std::cout << "]" << std::endl;
}

std::string GetMinimalErrorFilename(const FValuesMap &m, unsigned index) {
  double min_f = 1e8;
  std::string filename;
  for (FValuesMap::const_iterator it = m.begin();
       it != m.end(); ++it) {
    if (it->second[index] < min_f) {
      min_f = it->second[index];
      filename = it->first;
    }
  }
  // std::cout << filename << std::endl;
  return filename;
}
