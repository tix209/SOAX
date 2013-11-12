#ifndef SOAX_MULTISNAKE_H_
#define SOAX_MULTISNAKE_H_

#include "global.h"
#include "snake.h"
#include "junctions.h"


class QProgressBar;

namespace soax {

class SolverBank;

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

  void LoadParameters(const std::string &filename);
  void UpdateSnakeParameters();
  void SaveParameters(const std::string &filename) const;
  void PrintParameters() const;

  double intensity_scaling() const  {return intensity_scaling_;}
  void set_intensity_scaling(double scale) {intensity_scaling_ = scale;}
  bool intensity_scaled() const {return intensity_scaled_;}

  double sigma() const {return sigma_;}
  void set_sigma(double sigma) {sigma_ = sigma;}

  double ridge_threshold() const {return ridge_threshold_;}
  void set_ridge_threshold(double threshold) {ridge_threshold_ = threshold;}

  double foreground() const {return foreground_;}
  void set_foreground(double foreground) {foreground_ = foreground;}

  double background() const {return background_;}
  void set_background(double background) {background_ = background;}

  double stretch_factor() const {return stretch_factor_;}
  void set_stretch_factor(double factor) {stretch_factor_ = factor;}

  /*
   * Multiply the image intensity by a constant factor.
   */
  void ScaleImageIntensity();

  /*
   * Compute image gradient field for both snake initialization and
   * evolution.
   */
  void ComputeImageGradient();

  void InitializeSnakes();

  void SortSnakesOnLength(SnakeContainer &snakes);

  SnakeContainer &initial_snakes() {return initial_snakes_;}
  SnakeContainer &converged_snakes() {return converged_snakes_;}
  SnakeContainer &comparing_snakes1() {return comparing_snakes1_;}
  SnakeContainer &comparing_snakes2() {return comparing_snakes2_;}

  void SaveSnakes(const SnakeContainer &snakes,
                  const std::string &filename) const;

  void DeformSnakes(QProgressBar * progress_bar = NULL);
  void CutSnakesAtTJunctions();
  void GroupSnakes();


 private:
  typedef itk::Vector<bool, kDimension> BoolVectorType;
  typedef itk::Image<BoolVectorType, kDimension> BoolVectorImageType;
  typedef std::set<Snake *> SnakeSet;

  void SaveParameters(std::ofstream &outfile) const;

  /*
   * Initialize a new bool vector image used for snake initialization.
   */
  BoolVectorImageType::Pointer InitializeBoolVectorImage();

  /*
   * Scans the image gradient field to locate ridge points which is
   * significant (controlled by ridge_threshold_) and labeled in the
   * ridge_image.
   */
  void ScanGradient(BoolVectorImageType::Pointer ridge_image);

  /*
   * Generate candidate snake points.
   */
  void GenerateCandidates(BoolVectorImageType::Pointer ridge_image,
                          BoolVectorImageType::Pointer candidate_image,
                          unsigned direction);
  /*
   * Link candidate snake points into snakes.
   */
  void LinkCandidates(BoolVectorImageType::Pointer candidate_image,
                      unsigned direction);

  void LinkFromIndex(BoolVectorImageType::Pointer candidate_image,
                     BoolVectorImageType::IndexType &index,
                     unsigned direction);

  bool FindNextCandidate(BoolVectorImageType::Pointer candidate_image,
                         BoolVectorImageType::IndexType &index,
                         unsigned direction);


  static bool IsShorter(Snake *s1, Snake *s2) {
    return s1->length() < s2->length();
  }

  void AssignParameters(const std::string &name,
                        const std::string &value);

  double String2Double(const std::string &s);
  unsigned String2Unsigned(const std::string &s);

  void CutSnakes(SnakeContainer &seg);
  void ClearSnakeContainer(SnakeContainer &snakes);

  void LinkSegments(SnakeContainer &seg);
  void LinkFromSegment(Snake *s, SnakeContainer &seg,
                       SnakeSet &log, PointContainer &pc, bool &is_open);
  void LinkFromSegmentTip(SnakeTip *neighbor, PointContainer &pc,
                          bool &is_open, SnakeContainer &seg,
                          SnakeSet &log, bool from_head);
  void AddToPointContainer(PointContainer &pc, Snake *s,
                           bool is_head, bool from_head);
  void UpdateJunctions();
  unsigned GetNumberOfSnakesCloseToPoint(const PointType &p);

  std::string image_filename_;
  ImageType::Pointer image_;
  VectorImageType::Pointer external_force_;

  InterpolatorType::Pointer interpolator_;
  VectorInterpolatorType::Pointer vector_interpolator_;
  TransformType::Pointer transform_;
  SolverBank *solver_bank_;

  SnakeContainer initial_snakes_;
  SnakeContainer converged_snakes_;
  SnakeContainer comparing_snakes1_;
  SnakeContainer comparing_snakes2_;

  Junctions junctions_;




  /*
   * Intensity scale factor to normalize the intensity roughly inside
   * the range of [0, 1].
   */
  double intensity_scaling_;

  /*
   * True if image intensity is scaled.
   */
  bool intensity_scaled_;

  /*
   * Gaussian kernel size for smoothing the image.
   */
  double sigma_;

  /*
   * Ridge threshold for snake initialization. The bigger this value
   * is; the less snakes initialized.
   */
  double ridge_threshold_;

  /*
   * These two values specify the intensity range in which the snakes
   * are initialized. The stretching force is set to zero when the
   * snake tip intensity is below background_. When compute the local
   * stretch, voxels with intensity below background_ are not used as
   * samples.
   */
  double foreground_;
  double background_;

  /*
   * The distance between adjacent snake points the snakes try to maintain.
   */
  double desired_spacing_;

  /*
   * True if initialize snakes along z axis direction.
   */
  bool initialize_z_;

  /*
   * Minimum length for the final snake.
   */
  double minimum_length_;

  /*
   * Maximum number iterations allowed for a snake.
   */
  unsigned max_iterations_;

  /*
   * Change threshold for determining snake convergence. If every
   * snake point move a distance less than this threshold, then the
   * snake is converged.
   */
  double change_threshold_;

  /*
   * Period (# of iterations) of checking convergence.
   */
  unsigned check_period_;

  /*
   * Number of iterations performed on each click in step evolution mode.
   */
  unsigned iterations_per_press_;

  /*
   * Weight for first order continuity of snakes.
   */
  double alpha_;

  /*
   * Weight for second order continuity of snakes.
   */
  double beta_;

  /*
   * Step size of iteration.
   */
  double gamma_;

  /*
   * Weight for external forces.
   */
  double external_factor_;

  /*
   * Weight for stretch forces. It controls how much snakes stretch.
   */
  double stretch_factor_;

  /*
   * These thress parameters determine the sampling locations of local
   * shell in estimating the local background intensity near snake
   * tips.
   */
  int number_of_sectors_;
  int radial_near_;
  int radial_far_;

  /*
   * Number of points apart to compute the tangent vector at snake
   * tips.
   */
  unsigned delta_;

  /*
   * Distance threshold for determining overlap with other snakes.
   */
  double overlap_threshold_;

  /*
   * Distance threshold for determining junctions.
   */
  double grouping_distance_threshold_;

  /*
   * Number of points apart to compute the snake branch directions for
   * grouping.
   */
  unsigned grouping_delta_;

  /*
   * Angle threshold (in radians) for determining if two snake
   * branches are smooth enough to be linked together during grouping
   * process.
   */
  double direction_threshold_;

  /*
   * Flag of damping of stretching along z direction. If it is true,
   * the stretching forces are reduced if the tangent directions of
   * tips are along z direction.
   */
  bool damp_z_;


  DISALLOW_COPY_AND_ASSIGN(Multisnake);
};

} // namespace soax

#endif // SOAX_MULTISNAKE_H_
