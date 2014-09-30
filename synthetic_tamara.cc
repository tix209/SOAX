/** For Tamara synthetic image generation. PSF = (1, 1, 2). Gaussian noise sigma = 1.0.
 *  The image size is 150x150x150, and the offset for the coordinates is +60.0.
 */

#include <iostream>
#include "multisnake.h"

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: ./syn2 <curve-file> <output-image-file>" << std::endl;
    return -1;
  }

  soax::Multisnake ms;
  ms.LoadCurves(argv[1]);
  std::cout << ms.GetNumberOfComparingSnakes1() << " curves loaded." << std::endl;
  ms.PrintSnakes(ms.comparing_snakes1());
  ms.GenerateSyntheticTamara(argv[2]);
  return 0;
}
