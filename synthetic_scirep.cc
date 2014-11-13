/** Generate synthetic actin network images with shot noise. The image
 * SNR is controlled by the scaling factor S.
 */

#include <iostream>
#include <random>
#include "multisnake.h"


int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: ./syn3 <image-path> <snake-path> <output-dir>" << std::endl;
    return -1;
  }

  soax::Multisnake multisnake;
  multisnake.LoadImage(argv[1]);
  multisnake.LoadGroundTruthSnakes(argv[2]);
  std::cout << multisnake.GetNumberOfComparingSnakes1()
            << " ground truth snakes loaded." << std::endl;

  const unsigned offset = 0;
  const unsigned background = 30;
  const unsigned foreground = 90;
  const double scaling = 0.1;

  for (int i = 1; i <= 10; i++) {
    std::ostringstream buffer;
    buffer << "fg" << foreground << "-bg" << background << "-s" << scaling * i << ".tif";
    std::string output_file_path = std::string(argv[3]) + buffer.str();
    multisnake.GenerateSyntheticImageShotNoise(foreground, background, offset, scaling * i,
                                               output_file_path);
    // const std::string label("-diff");
    // std::string diff_filename = output_file_path;
    // diff_filename.insert(diff_filename.begin() + diff_filename.find_last_of('.') - 1,
    //                      label.begin(), label.end());
    // std::cout << diff_filename << std::endl;
    std::cout << output_file_path << std::endl;
  }
  return 0;
}
