/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file defines the SOAC (snake) class for SOAX.
 */


#ifndef SNAKE_H_
#define SNAKE_H_

#include <utility>
#include <vector>
#include <set>
#include "./global.h"

namespace soax {

class SolverBank;


class Snake {
 public:
  Snake(const PointContainer &points, bool is_open = true,
        bool is_grouping = false, ImageType::Pointer image = NULL,
        VectorImageType::Pointer external_force = NULL,
        InterpolatorType::Pointer interpolator = NULL,
        VectorInterpolatorType::Pointer vector_interpolator = NULL,
        TransformType::Pointer transform = NULL);

  bool open() const {return open_;}
  bool viable() const {return viable_;}
  double length() const {return length_;}
  double spacing() const {return spacing_;}
  bool converged() const {return converged_;}

  bool initial_state() const {return initial_state_;}
  void set_initial_state(bool initial) {initial_state_ = initial;}

  unsigned GetSize() const {return vertices_.size();}
  double GetX(unsigned i) const {return vertices_.at(i)[0];}
  double GetY(unsigned i) const {return vertices_.at(i)[1];}
  double GetZ(unsigned i) const {return vertices_.at(i)[2];}
  const PointType &GetPoint(unsigned i) const {return vertices_.at(i);}
  const PointType &GetHead() const {return vertices_.front();}
  const PointType &GetTail() const {return vertices_.back();}
  const PointType &GetTip(bool is_head) const;
  const PointContainer & vertices() const {return vertices_;}

  /*
   * Resample snake points to make them equally spaced close to the
   * specified spacing parameter (desired_spacing_). It will update
   * length_, spacing_ and viable_.
   */
  void Resample();

  void PrintSelf() const;
  void PrintVectorContainer(const VectorContainer &vc);

  void Trim(unsigned start, unsigned end);
  void ExtendHead(const PointType &p);
  void ExtendTail(const PointType &p);
  void TrimAndInsert(unsigned start, unsigned end, const PointType &p);

  // static void set_solver_bank(SolverBank *bank) {solver_bank_ = bank;}

  static double intensity_scaling() {return intensity_scaling_;}
  static void set_intensity_scaling(double scale) {
    intensity_scaling_ = scale;
  }

  static unsigned foreground() {return foreground_;}
  static void set_foreground(unsigned foreground) {
    foreground_ = foreground;
  }

  static unsigned background() {return background_;}
  static void set_background(unsigned background) {
    background_ = background;
  }

  static double desired_spacing() {return desired_spacing_;}
  static void set_desired_spacing(double spacing) {
    desired_spacing_ = spacing;
  }

  static double minimum_length() {return minimum_length_;}
  static void set_minimum_length(double length) {minimum_length_ = length;}

  static unsigned max_iterations() {return max_iterations_;}
  static void set_max_iterations(unsigned n) {max_iterations_ = n;}

  static double change_threshold() {return change_threshold_;}
  static void set_change_threshold(double v) {change_threshold_ = v;}

  static unsigned check_period() {return check_period_;}
  static void set_check_period(unsigned v) {check_period_ = v;}

  static unsigned iterations_per_press() {return iterations_per_press_;}
  static void set_iterations_per_press(unsigned n) {
    iterations_per_press_ = n;
  }

  // static void set_gamma(double g) {gamma_ = g;}

  static double external_factor() {return external_factor_;}
  static void set_external_factor(double f) {external_factor_ = f;}

  static double stretch_factor() {return stretch_factor_;}
  static void set_stretch_factor(double f) {stretch_factor_ = f;}

  static int number_of_sectors() {return number_of_sectors_;}
  static void set_number_of_sectors(int nsectors) {
    number_of_sectors_ = nsectors;
  }

  static int radial_near() {return radial_near_;}
  /*
   * Note radial_near_ cannot be less than 1 because the snake point
   * has to be included in the local foreground neighborhood.
   */
  static void set_radial_near(int rnear) {
    radial_near_ = rnear > 0 ? rnear : 1;
  }

  static int radial_far() {return radial_far_;}
  /*
   * Note radial_far must be greater than radial_near_.
   */
  static void set_radial_far(int rfar) {
    radial_far_ = rfar > 0 ? rfar : 2;
  }
  
  static int radial_save_foreground() {return radial_save_foreground_;}
  /*
   * Note radial_save_foreground_ cannot be less than 0
   */
  static void set_radial_save_foreground(int rforeground) {
    radial_save_foreground_ = rforeground >= 0 ? rforeground : 0;
  }

  static double z_spacing() {return z_spacing_;}
  static void set_z_spacing(double spacing) {
    z_spacing_ = spacing;
  }

  static unsigned delta() {return delta_;}
  static void set_delta(unsigned n) {delta_ = n;}

  static double overlap_threshold() {return overlap_threshold_;}
  static void set_overlap_threshold(double t) {
    overlap_threshold_ = t;
  }

  static double grouping_distance_threshold() {
    return grouping_distance_threshold_;
  }
  static void set_grouping_distance_threshold(double threshold) {
    grouping_distance_threshold_ = threshold;
  }

  static unsigned grouping_delta() {return grouping_delta_;}
  static void set_grouping_delta(unsigned d) {grouping_delta_ = d;}

  static double direction_threshold() {return direction_threshold_;}
  static void set_direction_threshold(double t) {direction_threshold_ = t;}

  static bool damp_z() {return damp_z_;}
  static void set_damp_z(bool d) {damp_z_ = d;}

  const SnakeContainer &subsnakes() const {return subsnakes_;}

  void Evolve(SolverBank *solver, const SnakeContainer &converged_snakes,
              unsigned max_iter, unsigned dim, const std::vector<std::vector<IndexPairContainer > >  &converged_snakes_grid);
  void EvolveWithTipFixed(SolverBank *solver, unsigned max_iter, unsigned dim);
  void UpdateHookedIndices();
  void CopySubSnakes(SnakeContainer &c);
  bool PassThrough(const PointType &p, double threshold) const;

  bool ComputeLocalSNRAtIndex(unsigned index, int radial_near, int radial_far,
                              double &local_snr) const;

  /*
   * Return the average intensity along the snake.
   */
  double ComputeIntensity() const;

  /*
   * Compute the mean background intensity around the tips. The sample points
   * are on a orthogonal plane at tips. The sampling region is a annulus
   * defined by "radial_near_" and "radial_far_".
   */
  double ComputeBackgroundMeanIntensity(unsigned index) const;
  double ComputeBackgroundMeanIntensity2d(unsigned index) const;

  double ComputeForegroundMeanIntensity(unsigned index) const;
  double ComputeForegroundMeanIntensity2d(unsigned index) const;

  InterpolatorOutputType InterpolateImageIntensity(PointType coords) const;

  /**
   * Return the average local SNR of the snake.
   */
  double ComputeSNR() const;

  bool IsInsideImage(const PointType &point,
                     unsigned dim = kDimension,
                     double padding = kBoundary) const;
 private:
  typedef std::vector<std::pair<double, double> > PairContainer;
  typedef std::set<unsigned> IndexSet;

  /*
   * Update length_ and compute length cumulative sum, stored in sums.
   */
  void UpdateLength(PairContainer *sums);

  /*
   * Compute new size for the current resample.
   */
  double ComputeNewSize(double spacing) const;

  /*
   * Update vertices_ by linearly interpolating the vertices based on
   * length cumulative sum.
   */
  void InterpolateVertices(const PairContainer *sums, unsigned new_size);

  bool IsConverged();
  void CheckSelfIntersection();
  bool TipsStopAtSameLocation();
  void TryInitializeFromPart(PointIterator it1, PointIterator it2,
                             bool is_open);

  void HandleHeadOverlap(const SnakeContainer &converged_snakes, const std::vector<std::vector<IndexPairContainer > > &converged_snakes_grid);
  void HandleTailOverlap(const SnakeContainer &converged_snakes, const std::vector<std::vector<IndexPairContainer > > &converged_snakes_grid);

  bool HeadIsFixed() {return fixed_head_[0] > 0;}
  bool TailIsFixed() {return fixed_tail_[0] > 0;}

  PointIterator CheckHeadOverlap(PointIterator const &start,
                                 const SnakeContainer &converged_snakes, const std::vector<std::vector<IndexPairContainer > >  &converged_snakes_grid);
  PointIterator CheckTailOverlap(PointIterator const &start,
                                 const SnakeContainer &converged_snakes, const std::vector<std::vector<IndexPairContainer > >  &converged_snakes_grid);

  bool VertexOverlap(const PointType &p,
                     const SnakeContainer &converged_snakes, const std::vector<std::vector<IndexPairContainer > >  &converged_snakes_grid);

  void FindHookedSnakeAndIndex(const PointType &p,
                               const SnakeContainer &converged_snakes,
                               Snake * &s, unsigned &index, const std::vector<std::vector<IndexPairContainer > >  &converged_snakes_grid);
  double FindClosestIndexTo(const PointType &p, unsigned &ind);

  void IterateOnce(SolverBank *solver, unsigned dim);

  void ComputeRHSVector(double gamma, VectorContainer &rhs, unsigned dim);
  void AddExternalForce(VectorContainer &rhs, unsigned dim);
  void AddStretchingForce(VectorContainer &rhs, unsigned dim);
  void AddVerticesInfo(double gamma, VectorContainer &rhs);

  void UpdateHeadTangent();
  void UpdateTailTangent();

  double ComputeLocalStretch(unsigned index, unsigned dim);

  double ComputePodX(double x, const VectorType &tvec,
                     double dist, bool plus_root) const;
  double ComputePodY(double y, const VectorType &tvec,
                     double dist, bool plus_root) const;
  bool CheckOrthogonality(const VectorType &vec1,
                          const VectorType &vec2) const;

  void GetStartingRadialDirection(VectorType &direction,
                                  const VectorType &normal,
                                  const PointType &vertex) const;
  unsigned GetPrincipalIndex(const VectorType &vec) const;
  void ComputeSamplePoint(PointType &point,
                          const PointType &origin,
                          const VectorType &radial,
                          const VectorType &normal,
                          int d, int s) const;

  void CheckBodyOverlap(const SnakeContainer &converged_snakes, const std::vector<std::vector<IndexPairContainer > >  &converged_snakes_grid);

  void AddJunctionIndex(unsigned index);

  double ComputeLocalForegroundMean(unsigned index, int radial_near) const;
  bool ComputeLocalBackgroundMeanStd(unsigned index,
                                     int radial_near, int radial_far,
                                     double &mean, double &std) const;
  VectorType ComputeUnitTangentVector(unsigned index) const;



  PointContainer vertices_;
  bool open_;
  bool grouping_;
  ImageType::Pointer image_;
  VectorImageType::Pointer external_force_;
  InterpolatorType::Pointer interpolator_;
  VectorInterpolatorType::Pointer vector_interpolator_;
  TransformType::Pointer transform_;

  PointContainer last_vertices_;
  bool viable_;

  /*
   * The snake spacing is 1/(kMimimumEvolvingSize-1) = 0.25 if
   * initial_state_ is True; otherwise desired_spacing_ is used.
   */
  bool initial_state_;
  bool converged_;
  bool final_;

  double length_;
  double spacing_;

  /*
   * The average intensity along this snake.
   */
  double intensity_;
  unsigned iterations_;

  VectorType head_tangent_;
  VectorType tail_tangent_;

  PointType fixed_head_;
  PointType fixed_tail_;

  Snake *head_hooked_snake_;
  Snake *tail_hooked_snake_;
  unsigned head_hooked_index_;
  unsigned tail_hooked_index_;

  /*
   * Snake can be devided into several subsnakes during its evolution.
   */
  SnakeContainer subsnakes_;
  IndexSet junction_indices_;

  // static SolverBank *solver_bank_;

  static double intensity_scaling_;
  static unsigned foreground_;
  static unsigned background_;

  static double desired_spacing_;

  /*
   * Maximum number iterations allowed for a snake.
   */
  static unsigned max_iterations_;

  /*
   * Change threshold for determining snake convergence. If every
   * snake point move a distance less than this threshold, then the
   * snake is converged.
   */
  static double change_threshold_;

  /*
   * Period (# of iterations) of checking convergence.
   */
  static unsigned check_period_;

  static unsigned iterations_per_press_;
  /*
   * Minimum length for the final snake.
   */
  static double minimum_length_;

  // static double gamma_;
  /*
   * Weight for external forces.
   */
  static double external_factor_;

  /*
   * Weight for stretch forces.
   */
  static double stretch_factor_;

  /*
   * These thress parameters determine the sampling locations of local
   * shell in estimating the local background intensity near snake
   * tips.
   */
  static int number_of_sectors_;
  static int radial_near_;
  static int radial_far_;
  static int radial_save_foreground_;

  /*
   * Number of points apart to compute the tangent vector at snake
   * tips.
   */
  static unsigned delta_;

  /*
   * Distance threshold for determining overlap with other snakes.
   */
  static double overlap_threshold_;

  /*
   * Distance threshold for determining junctions.
   */
  static double grouping_distance_threshold_;

  /*
   * Number of points apart to compute the snake branch directions for
   * grouping.
   */
  static unsigned grouping_delta_;

  /*
   * Angle threshold (in radians) for determining if two snake
   * branches are smooth enough to be linked together during grouping
   * process.
   */
  static double direction_threshold_;

  /*
   * Flag of damping of stretching along z direction. If it is true,
   * the stretching forces are reduced if the tangent directions of
   * tips are along z direction.
   */
  static bool damp_z_;
  /*
   * Image boundary size in pixels.
   */
  static const double kBoundary;

  /**
   * Voxel size relative to x/y.
   */
  static double z_spacing_;

  DISALLOW_COPY_AND_ASSIGN(Snake);
};

}  // namespace soax

#endif  // SNAKE_H_
