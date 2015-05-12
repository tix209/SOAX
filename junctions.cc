/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file implements the network junctions class for SOAX.
 */

#include <cassert>
#include <fstream>
#include <algorithm>
#include "./junctions.h"
#include "./snake.h"

namespace soax {
Junctions::Junctions() {}

Junctions::~Junctions() {
  this->Reset();
}


void Junctions::DeleteTipSetContainer() {
  if (tip_sets_.empty()) return;

  for (TipSetContainer::iterator it = tip_sets_.begin();
       it != tip_sets_.end(); ++it) {
    delete *it;
  }
  tip_sets_.clear();
}

void Junctions::DeleteTipContainer() {
  if (tips_.empty()) return;
  for (TipContainer::iterator it = tips_.begin(); it != tips_.end(); ++it) {
    delete *it;
  }
  tips_.clear();
}

void Junctions::Reset() {
  this->DeleteTipSetContainer();
  this->DeleteTipContainer();
  junction_points_.clear();
}

void Junctions::Initialize(const SnakeContainer &seg) {
  for (SnakeConstIterator it = seg.begin(); it != seg.end(); ++it) {
    SnakeTip * head = new SnakeTip(*it, true);
    SnakeTip * tail = new SnakeTip(*it, false);

    tips_.push_back(head);
    tips_.push_back(tail);
  }
}

void Junctions::Union() {
  assert(!tips_.empty());

  for (TipContainer::iterator it = tips_.begin();
       it != tips_.end(); ++it) {
    this->AddToTipSets(*it);
  }
  // Get rid of both tips included in a snaketipset
  TipContainer singleton_tips;
  for (TipSetContainer::iterator it = tip_sets_.begin();
       it != tip_sets_.end(); ++it) {
    (*it)->Purge(singleton_tips);
  }
  for (TipContainer::iterator it = singleton_tips.begin();
       it != singleton_tips.end(); ++it) {
    this->AddNewTipSet(*it);
  }

  // Get the filament junction locations
  for (TipSetContainer::iterator it = tip_sets_.begin();
       it != tip_sets_.end(); ++it) {
    (*it)->UpdatePosition();
    if ((*it)->GetDegree() > 1)
      junction_points_.push_back((*it)->position());
  }
}

void Junctions::AddToTipSets(SnakeTip *t) {
  TipSetContainer close_tip_sets;
  // Find any TipSets that is close to t and put them to close_tip_sets
  this->FindCloseTipSets(t, close_tip_sets);
  // Get consolidated TipSet by merging TipSets which are all close to t
  SnakeTipSet *ts = this->MergeCloseTipSets(close_tip_sets);
  // Add t to the consolidated TipSet
  this->JoinTipSet(t, ts);
}

void Junctions::FindCloseTipSets(SnakeTip *t,
                                 TipSetContainer &close_tip_sets) const {
  if (tip_sets_.empty()) return;

  for (TipSetContainer::const_iterator it = tip_sets_.begin();
       it != tip_sets_.end(); ++it) {
    double dist = (*it)->ComputeDistance(t);
    if (dist < Snake::grouping_distance_threshold()) {
      close_tip_sets.push_back(*it);
    }
  }
}

SnakeTipSet * Junctions::MergeCloseTipSets(
    const TipSetContainer &close_tip_sets) {
  if (close_tip_sets.empty()) {
    return NULL;
  } else if (close_tip_sets.size() == 1) {
    return close_tip_sets.front();
  } else {
    SnakeTipSet *ts = close_tip_sets.front();
    for (TipSetContainer::const_iterator it = close_tip_sets.begin()+1;
         it != close_tip_sets.end(); ++it) {
      ts->Combine(*it);
      this->RemoveTipSet(*it);
    }

    return ts;
  }
}


void Junctions::RemoveTipSet(SnakeTipSet *ts) {
  TipSetContainer::iterator it = std::find(tip_sets_.begin(),
                                           tip_sets_.end(),
                                           ts);
  if (it != tip_sets_.end()) {
    tip_sets_.erase(it);
  } else {
    std::cerr << "Junctions::RemoveTipSet: ts not found!" << std::endl;
  }
}


void Junctions::JoinTipSet(SnakeTip *t, SnakeTipSet *ts) {
  if (!ts)
    this->AddNewTipSet(t);
  else
    ts->Add(t);
}


void Junctions::AddNewTipSet(SnakeTip *t) {
  SnakeTipSet *ts = new SnakeTipSet(t);
  tip_sets_.push_back(ts);
}

double Junctions::FindTipSetToJoin(SnakeTip *t, SnakeTipSet * &ts) {
  double min_dist = tip_sets_.front()->ComputeDistance(t);
  ts = tip_sets_.front();
  TipSetContainer::iterator it = tip_sets_.begin();
  it++;
  while (it != tip_sets_.end()) {
    double d = (*it)->ComputeDistance(t);
    if (d < min_dist) {
      min_dist = d;
      ts = *it;
    }
    ++it;
  }

  return min_dist;
}


void Junctions::Configure() {
  assert(!tip_sets_.empty());
  for (TipSetContainer::iterator it = tip_sets_.begin();
       it != tip_sets_.end(); ++it) {
    (*it)->Configure();
  }
}


void Junctions::PrintTips() const {
  for (TipContainer::const_iterator it = tips_.begin();
       it != tips_.end(); ++it) {
    (*it)->PrintSelf();
  }
}

void Junctions::PrintTipSets() const {
  for (TipSetContainer::const_iterator it = tip_sets_.begin();
       it != tip_sets_.end(); ++it) {
    (*it)->PrintSelf();
  }
}

SnakeTip * Junctions::FindSnakeTip(Snake *s, bool is_head) const {
  for (TipContainer::const_iterator it = tips_.begin();
       it != tips_.end(); ++it) {
    if ((*it)->snake() == s && (*it)->is_head() == is_head)
      return *it;
  }
  return NULL;
}

void Junctions::PrintJunctionPoints(const std::string &filename) const {
  std::ofstream outfile;
  outfile.open(filename.c_str(), std::ios::out | std::ios::app);
  if (!outfile.is_open()) {
    std::cerr << "Junctions::PrintJunctionPoints: Couldn't open file: "
              << outfile << std::endl;
    return;
  }

  for (PointContainer::const_iterator it = junction_points_.begin();
       it != junction_points_.end(); ++it) {
    outfile << *it << std::endl;
  }
  outfile.close();
}


void Junctions::RemoveJunction(const PointType &p) {
  assert(!junction_points_.empty());
  junction_points_.erase(std::find(junction_points_.begin(),
                                   junction_points_.end(),
                                   p));
}

void Junctions::RemoveJunctions(const std::vector<PointType> &pts) {
  for (int i = 0; i < pts.size(); i++) {
    RemoveJunction(pts[i]);
  }
}

void Junctions::ClearJunctionPoints() {
  junction_points_.clear();
}

double Junctions::ComputeMinDistance(const PointType &p,
                                     const PointContainer &pts) {
  double min_d = p.EuclideanDistanceTo(pts.front());
  for (unsigned i = 1; i < pts.size(); ++i) {
    double d = p.EuclideanDistanceTo(pts.at(i));
    if (d < min_d)
      min_d = d;
  }
  return min_d;
}

}  // namespace soax
