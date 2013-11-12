#ifndef SOAX_SNAKE_TIP_H_
#define SOAX_SNAKE_TIP_H_

#include "global.h"

namespace soax {

class SnakeTip {
 public:
  SnakeTip(Snake *s, bool is_head);

  PointType tip() const {return tip_;}
  Snake * snake() const {return snake_;}
  bool is_head() const {return is_head_;}
  SnakeTip * neighbor() const {return neighbor_;}
  void set_neighbor(SnakeTip *neighbor) {neighbor_ = neighbor;}

  double DistanceTo(SnakeTip *t) const;

  void Link(SnakeTip *t) {neighbor_ = t;}

  VectorType GetDirection() const;
  static double ComputeAngle(const SnakeTip *t1, const SnakeTip *t2);

  void PrintSelf() const;


 private:
  Snake *snake_;
  bool is_head_;
  PointType tip_;
  SnakeTip *neighbor_;

  DISALLOW_COPY_AND_ASSIGN(SnakeTip);
};
} // namespace soax
#endif //SOAX_SNAKE_TIP_H_
