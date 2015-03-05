/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file implements the SOAC tip set class for the network configuration
 * procedure in SOAX.
 */

#include <cassert>
#include "./snake_tip_set.h"
#include "./snake.h"

namespace soax {

SnakeTipSet::SnakeTipSet(SnakeTip *t) {
  this->Add(t);
}

void SnakeTipSet::Add(SnakeTip *t) {
  tips_.push_back(t);
}

void SnakeTipSet::Combine(SnakeTipSet *ts) {
  for (TipList::const_iterator it = ts->tips_.begin();
       it != ts->tips_.end(); ++it) {
    this->Add(*it);
  }
}

void SnakeTipSet::UpdatePosition() {
  degree_ = tips_.size();
  assert(!tips_.empty());
  for (unsigned i = 0; i < kDimension; ++i) {
    double sum = 0;
    for (TipList::iterator it = tips_.begin();
         it != tips_.end(); ++it) {
      sum += (*it)->snake()->GetTip((*it)->is_head())[i];
    }
    position_[i] = sum / degree_;
  }
}

double SnakeTipSet::ComputeDistance(SnakeTip *t) const {
  assert(!tips_.empty());
  double min_d = t->DistanceTo(tips_.front());
  TipList::const_iterator it = tips_.begin();
  it++;
  while (it != tips_.end()) {
    double d = t->DistanceTo(*it);
    if (d < min_d) {
      min_d = d;
    }

    it++;
  }
  return min_d;
}

void SnakeTipSet::Configure() {
  if (tips_.size() < 2) return;

  SnakeTip *t1 = NULL;
  SnakeTip *t2 = NULL;

  double angle = this->FindSmoothestPair(t1, t2);
  if (angle < Snake::direction_threshold()) {
    return;
  } else {
    t1->Link(t2);
    t2->Link(t1);
    tips_.remove(t1);
    tips_.remove(t2);
    this->Configure();
  }
}

double SnakeTipSet::FindSmoothestPair(SnakeTip * &t1, SnakeTip * &t2) {
  double max_angle = 0.0;
  for (TipList::iterator i = tips_.begin(); i != tips_.end(); ++i) {
    for (TipList::iterator j = i; j != tips_.end(); ++j) {
      if (j != i) {
        double angle = SnakeTip::ComputeAngle(*i, *j);
        if (angle > max_angle) {
          max_angle = angle;
          t1 = *i;
          t2 = *j;
        }
      }
    }
  }
  return max_angle;
}

void SnakeTipSet::PrintSelf() const {
  std::cout << "************* SnakeTipSet *************" << std::endl;
  std::cout << degree_ << "\t" << position_ << std::endl;
  std::cout << "Snake tips: " << std::endl;
  for (TipList::const_iterator it = tips_.begin();
       it != tips_.end(); ++it) {
    (*it)->PrintSelf();
  }
  std::cout << std::endl;
}

void SnakeTipSet::Purge(TipContainer &tips) {
  for (TipList::iterator i = tips_.begin(); i != tips_.end(); ++i) {
    for (TipList::iterator j = i; j != tips_.end(); ++j) {
      if (j != i) {
        if ((*i)->snake() == (*j)->snake()) {
          double dist_i = this->ComputeTotalDistanceToOthers(*i);
          double dist_j = this->ComputeTotalDistanceToOthers(*j);
          if (dist_i < dist_j) {
            tips.push_back(*j);
          } else {
            tips.push_back(*i);
          }
        }
      }
    }
  }
  for (TipContainer::const_iterator it = tips.begin();
       it != tips.end(); ++it) {
    tips_.remove(*it);
  }
}

double SnakeTipSet::ComputeTotalDistanceToOthers(SnakeTip *t) const {
  double total = 0.0;
  TipList::const_iterator it = tips_.begin();
  while (it != tips_.end()) {
    total += t->DistanceTo(*it);
    it++;
  }
  return total;
}

}  // namespace soax
