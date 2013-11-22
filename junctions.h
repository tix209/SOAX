#ifndef SOAX_JUNCTIONS_H_
#define SOAX_JUNCTIONS_H_

#include "global.h"
#include "snake_tip_set.h"

namespace soax {

class Junctions {
 public:
  Junctions();

  ~Junctions();

  const PointContainer &junction_points() const {return junction_points_;}
  void set_junction_points(const PointContainer &p) {junction_points_ = p;}
  void Initialize(const SnakeContainer &seg);
  void Union();
  void Configure();

  double ComputeDistance(SnakeTip *t);

  void PrintTips() const;
  void PrintTipSets() const;
  void PrintJunctionPoints(const std::string &filename) const;
  SnakeTip * FindSnakeTip(Snake *s, bool is_head) const;
  void Reset();

  /*
   * Remove junction specified by its location p. Note only
   * junction_points_ are changed.
   */
  void RemoveJunction(const PointType &p);

  void ClearJunctionPoints();

  /*
   * Remove junction points that are not in tips.
   */
  void UpdateJunctionPoints(const PointContainer &tips);


 private:
  // Junctions(const Junctions &j) {}
  // Junctions & operator=(const Junctions &j) {return *this;}

  typedef std::vector<SnakeTipSet *> TipSetContainer;
  typedef std::vector<SnakeTip *> TipContainer;

  void DeleteTipSetContainer();
  void DeleteTipContainer();
  void AddToTipSets(SnakeTip *t);

  void FindCloseTipSets(SnakeTip *t, TipSetContainer &close_tip_sets) const;

  SnakeTipSet * MergeCloseTipSets(const TipSetContainer &close_tip_sets);

  void RemoveTipSet(SnakeTipSet *ts);
  void JoinTipSet(SnakeTip *t, SnakeTipSet *ts);

  void AddNewTipSet(SnakeTip *t);


  /*
   * Returns the smallest distance from a snake tip set to t. ts is
   * set to the nearest snake tip set.
   */
  double FindTipSetToJoin(SnakeTip *t, SnakeTipSet * &ts);


  double ComputeMinDistance(const PointType &p, const PointContainer &pts);

  TipSetContainer tip_sets_;
  TipContainer tips_;
  PointContainer junction_points_;

  DISALLOW_COPY_AND_ASSIGN(Junctions);
};

} // namespace soax
#endif //SOAX_JUNCTIONS_H_
