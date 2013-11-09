#ifndef SOAX_MULTISNAKE_H_
#define SOAX_MULTISNAKE_H_

#include "global.h"

namespace soax {

/*
 * Multiple open snake class.
 */
class Multisnake {
 public:
  Multisnake();
  ~Multisnake();

  /*
   * Load the image and set image_filename_.
   */
  void LoadImage(const std::string &filename);
  ImageType::Pointer image() const {return image_;}

 private:
  std::string image_filename_;
  ImageType::Pointer image_;

  DISALLOW_COPY_AND_ASSIGN(Multisnake);
};

} // namespace soax

#endif // SOAX_MULTISNAKE_H_
