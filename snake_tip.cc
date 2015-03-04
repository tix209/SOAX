/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file implements the SOAC tip class for the network configuration
 * procedure in SOAX.
 */


#include <cmath>
#include "./snake_tip.h"
#include "./snake.h"

namespace soax {

SnakeTip::SnakeTip(Snake *s, bool is_head) :
    snake_(s), is_head_(is_head), neighbor_(NULL) {
  tip_ = s->GetTip(is_head);
}

double SnakeTip::DistanceTo(SnakeTip *t) const {
  if (t == this) return 0.0;
  return tip_.EuclideanDistanceTo(t->tip());
}

double SnakeTip::ComputeAngle(const SnakeTip *t1, const SnakeTip *t2) {
  VectorType v1 = t1->GetDirection();
  VectorType v2 = t2->GetDirection();

  return acos(v1 * v2);
}

VectorType SnakeTip::GetDirection() const {
  unsigned delta = snake_->GetSize() > Snake::grouping_delta() ?
      Snake::grouping_delta() : snake_->GetSize();
  VectorType v;
  if (is_head_) {
    v = snake_->GetPoint(delta/2) - snake_->GetPoint(delta-1);
    // v = snake_->GetHead() - snake_->GetPoint(delta-1);
  } else {
    v = snake_->GetPoint(snake_->GetSize()-delta/2) -
        snake_->GetPoint(snake_->GetSize()-delta);
    // v = snake_->GetTail() - snake_->GetPoint(snake_->GetSize()-delta);
  }
  v.Normalize();
  return v;
}


void SnakeTip::PrintSelf() const {
  std::cout << "\n~~~~~~ SnakeTip ~~~~~~" << std::endl;
  std::cout << this << std::endl;
  std::cout << "snake: " << snake_ << "\t" << is_head_ << std::endl;
  std::cout << tip_ << std::endl;
  std::cout << neighbor_ << std::endl;
}

}  // namespace soax
