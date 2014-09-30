/** For Tamara synthetic image generation. PSF = (1, 1, 2). Gaussian noise sigma = 1.0.
 *  The image size is 80x150x80.
 */

#include <iostream>
#include "multisnake.h"

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: ./syn2 <curve-file> <output-image-file>" << std::endl;
    return -1;
  }

  soax::Multisnake ms;
  double coordinates_offset[3] = {40, 70, 40};
  ms.LoadCurves(argv[1], coordinates_offset);
  std::cout << ms.GetNumberOfComparingSnakes1() << " curves loaded." << std::endl;
  // ms.PrintSnakes(ms.comparing_snakes1());
  double foreground = 40;
  double background = 230;
  ms.GenerateSyntheticTamara(foreground, background, argv[2]);
  return 0;
}
