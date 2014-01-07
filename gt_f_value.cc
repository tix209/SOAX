#include <iostream>
#include <cstdlib>
#include "multisnake.h"

/*
 * Compute the F-function value of a ground truth snake. The low SNR
 * threshold is t * image_snr, there the image_snr is estimated by
 * Otsu's method. The penalizing factor c is empirically chosen to 1.5.
 */

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: gtf <image_path> <gt_snake_path>" << std::endl;
    return EXIT_FAILURE;
  }
  std::string image_path = argv[1];
  std::string gt_snake_path = argv[2];

  soax::Multisnake multisnake;
  multisnake.LoadImage(image_path);
  double snr = multisnake.ComputeImageSNR();
  std::cout << "Image SNR: " << snr << std::endl;

  multisnake.LoadGroundTruthSnakes(gt_snake_path);
  std::cout << multisnake.GetNumberOfComparingSnakes1() << " snakes loaded."
            << std::endl;
  double threshold = 0.6 * snr;
  double penalizer = 1.5;
  double fvalue = multisnake.ComputeGroundTruthFValue(threshold, penalizer);
  std::cout << "F-value: " << fvalue << std::endl;
  return EXIT_SUCCESS;
}
