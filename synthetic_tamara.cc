/** For Tamara synthetic image generation. PSF = (1, 1, 2). Gaussian noise sigma = 1.0.
 *  The image size is 80x150x80.
 */

#include <iostream>
#include "multisnake.h"
#include "boost/filesystem.hpp"

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: ./syn2 <curve-dir> <output-image-dir>" << std::endl;
    return -1;
  }

  namespace fs = boost::filesystem;
  fs::path indir(argv[1]);
  // fs::path outdir(argv[2]);

  typedef std::vector<fs::path> Paths;
  Paths infiles;
  std::copy(fs::directory_iterator(indir), fs::directory_iterator(),
            back_inserter(infiles));
  std::sort(infiles.begin(), infiles.end());

  double coordinates_offset[3] = {40, 70, 40};
  const double foreground = 40;
  const double background = 230;

  for (Paths::const_iterator it(infiles.begin()); it != infiles.end(); ++it) {
    soax::Multisnake ms;
    std::cout << *it << std::endl;
    ms.LoadCurves(it->string().c_str(), coordinates_offset);
    std::cout << ms.GetNumberOfComparingSnakes1() << " curves loaded." << std::endl;

    std::string outfile = std::string(argv[2]) + it->filename().string();
    size_t dot_pos = outfile.find_last_of(".");
    outfile.replace(dot_pos + 1, 3, "mha");
    ms.GenerateSyntheticTamara(foreground, background, outfile.c_str());

    std::cout << outfile << std::endl;
  }

  return 0;
}
