#ifndef SOAX_DROPLET_INFO_H_
#define SOAX_DROPLET_INFO_H_

#include "global.h"

namespace soax {

class DropletInfo {
 public:
  DropletInfo();
  DropletInfo(double center_x, double center_y,
              double center_z, double radius, unsigned type);
  PointType center() const {return center_;}
  double radius() const {return radius_;}
  unsigned type() const {return type_;}

 private:
  PointType center_;
  double radius_;
  unsigned type_;
};

} // namespace soax
#endif // SOAX_DROPLET_INFO_H_
