#include "droplet_info.h"

namespace soax {

DropletInfo::DropletInfo() {}

DropletInfo::DropletInfo(double center_x, double center_y,
                         double center_z, double radius, unsigned type) {
  center_[0] = center_x;
  center_[1] = center_y;
  center_[2] = center_z;
  radius_ = radius;
  type = type_;
}
}
