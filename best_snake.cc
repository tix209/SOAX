#include <iostream>
#include <cstdlib>
#include <vector>
#include <map>
#include <fstream>
// #include <iomanip>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "multisnake.h"
// #include "utility.h"

typedef std::map<std::string, std::vector<double> > FValuesMap;
typedef std::map<double, std::string> TFilenameMap;
std::vector<double> ComputeFValueVector(soax::Multisnake &ms, int rnear, int rfar,
                                        double start, double end, double step);
void PrintFMap(const FValuesMap &fmap);
void PrintVector(const std::vector<double> &v);
void PrintTFilenameMap(const TFilenameMap &tfm, std::ostream &os);
std::string GetMinimalErrorFilename(const FValuesMap &m, unsigned index);



int main (int argc, char **argv) {
  if (argc < 6) {
    std::cerr << "./best_soacs <image_path> <snake_dir> "
        "<radial_near> <radial_far> <output_path> " << std::endl;
    return EXIT_FAILURE;
  }

  namespace fs = boost::filesystem;
  std::string image_path = argv[1];
  fs::path snake_dir(argv[2]);

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
      double t_end = 5.01;
      double t_step = 0.2;
      int rnear = std::atoi(argv[3]);
      int rfar = std::atoi(argv[4]);
      for (Paths::const_iterator it = sorted_snakes_path.begin();
           it != sorted_snakes_path.end(); ++it) {
        std::cout << it->filename() << std::endl;
        ms.LoadConvergedSnakes(it->string());
        fmap[it->filename().string()] = ComputeFValueVector(
            ms, rnear, rfar, t_start, t_end, t_step);
      }

      unsigned num_entries = static_cast<unsigned>((t_end - t_start) / t_step)+1;

      TFilenameMap t_filenames;
      for (unsigned i = 0; i < num_entries; i++) {
        std::string filename = GetMinimalErrorFilename(fmap, i);
        double t = t_start + t_step * i;
        t_filenames[t] = filename;
        // std::cout << filename << std::endl;
      }

      // PrintFMap(fmap);
      std::ofstream outfile(argv[5]);
      PrintTFilenameMap(t_filenames, outfile);
    } else {
      std::cout << snake_dir << " does not exist." << std::endl;
    }
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}


std::vector<double> ComputeFValueVector(soax::Multisnake &ms, int rnear, int rfar,
                                        double start, double end, double step) {
  const double c = 1.1;
  // const int rnear = 4;
  // const int rfar = 12;
  std::vector<double> v;
  for (double t = start; t < end; t += step) {
    soax::DataContainer snrs;
    ms.ComputeResultSnakesLocalSNRs(rnear, rfar, snrs);
    v.push_back(ms.ComputeFValue(snrs, t, c));
  }
  return v;
}

void PrintTFilenameMap(const TFilenameMap &tfm, std::ostream &os) {
  for (TFilenameMap::const_iterator it = tfm.begin();
       it != tfm.end(); ++it) {
    os << it->first << ": " << it->second << std::endl;
  }
}

void PrintFMap(const FValuesMap &fmap) {
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
