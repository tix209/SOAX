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
  void Reset();
  /*
   * Load the image and set image_filename_.
   */
  void LoadImage(const std::string &filename);

  /*
   * Resample and save as an isotropic 16-bit image.
   */
  void SaveAsIsotropicImage(const std::string &filename, double z_spacing);

  ImageType::Pointer image() const {return image_;}

  void LoadParameters(const std::string &filename);
  void UpdateSnakeParameters();
  void SaveParameters(const std::string &filename) const;
  void WriteParameters(std::ostream &os) const;

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

  bool initialize_z() const {return initialize_z_;}
  void set_initialize_z(bool init_z) {initialize_z_ = init_z;}

  SolverBank *solver_bank() const {return solver_bank_;}

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

  unsigned GetNumberOfInitialSnakes() const {return initial_snakes_.size();}
  unsigned GetNumberOfConvergedSnakes() const {
    return converged_snakes_.size();
  }
  unsigned GetNumberOfComparingSnakes1() const {
    return comparing_snakes1_.size();
  }
  unsigned GetNumberOfComparingSnakes2() const {
    return comparing_snakes2_.size();
  }

  const SnakeContainer &initial_snakes() const {return initial_snakes_;}
  const SnakeContainer &converged_snakes() const {return converged_snakes_;}
  const SnakeContainer &comparing_snakes1() const {return comparing_snakes1_;}
  const SnakeContainer &comparing_snakes2() const {return comparing_snakes2_;}

  void SaveConvergedSnakesAsJFilamentFormat(
      const std::string &filename) const {
    this->SaveJFilamentSnakes(converged_snakes_, filename);
  }

  void DeformSnakes(QProgressBar * progress_bar = NULL);
  void CutSnakesAtTJunctions();
  void GroupSnakes();

  const PointContainer &GetJunctions() const {
    return junctions_.junction_points();
  }

  void LoadConvergedSnakes(const std::string &filename) {
    this->LoadSnakes(filename, converged_snakes_);
  }

  void LoadGroundTruthSnakes(const std::string &filename) {
    this->LoadJFilamentSnakes(filename, comparing_snakes1_);
  }

  void LoadComparingSnakes1(const std::string &filename) {
    this->LoadSnakes(filename, comparing_snakes1_);
  }

  void LoadComparingSnakes2(const std::string &filename) {
    this->LoadSnakes(filename, comparing_snakes2_);
  }

  void PrintSnakes(const SnakeContainer &snakes) const;

  void SaveSnakes(const SnakeContainer &snakes,
                  const std::string &filename) const;

  void EvaluateByVertexErrorHausdorffDistance(
      const std::string &snake_path, const std::string &filename) const;
  void EvaluateByFFunction(double threshold, double penalizer,
                           int radial_near, int radial_far,
                           const std::string &snake_path,
                           const std::string &filename) const;

  void PrintGroundTruthLocalSNRValues(int radial_near, int radial_far);

 private:
  typedef itk::Vector<bool, kDimension> BoolVectorType;
  typedef itk::Image<BoolVectorType, kDimension> BoolVectorImageType;
  typedef std::set<Snake *> SnakeSet;


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

  void LoadSnakes(const std::string &filename, SnakeContainer &snakes);
  void LoadJFilamentSnakes(const std::string &filename,
                           SnakeContainer &snakes);
  void LoadPoint(const std::string &s, PointContainer &c);
  void SaveJFilamentSnakes(const SnakeContainer &snakes, const std::string &filename) const;

  void ComputeErrorFromSnakesToComparingSnakes(DataContainer &errors) const;
  void ComputeErrorFromComparingSnakesToSnakes(DataContainer &errors) const;
  double ComputeShortestDistance(const PointType &p,
                                 const SnakeContainer &snakes) const;


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

  // /*
  //  * The distance between adjacent snake points the snakes try to maintain.
  //  */
  // double desired_spacing_;

  /*
   * True if initialize snakes along z axis direction.
   */
  bool initialize_z_;

  // /*
  //  * Minimum length for the final snake.
  //  */
  // double minimum_length_;

  // /*
  //  * Maximum number iterations allowed for a snake.
  //  */
  // unsigned max_iterations_;

  // /*
  //  * Change threshold for determining snake convergence. If every
  //  * snake point move a distance less than this threshold, then the
  //  * snake is converged.
  //  */
  // double change_threshold_;

  // /*
  //  * Period (# of iterations) of checking convergence.
  //  */
  // unsigned check_period_;

  /*
   * Number of iterations performed on each click in step evolution mode.
   */
  // unsigned iterations_per_press_;

  // /*
  //  * Weight for first order continuity of snakes.
  //  */
  // double alpha_;

  // /*
  //  * Weight for second order continuity of snakes.
  //  */
  // double beta_;

  // /*
  //  * Step size of iteration.
  //  */
  // double gamma_;

  // /*
  //  * Weight for external forces.
  //  */
  // double external_factor_;

  // /*
  //  * Weight for stretch forces. It controls how much snakes stretch.
  //  */
  // double stretch_factor_;

  // /*
  //  * These thress parameters determine the sampling locations of local
  //  * shell in estimating the local background intensity near snake
  //  * tips.
  //  */
  // int number_of_sectors_;
  // int radial_near_;
  // int radial_far_;

  // /*
  //  * Number of points apart to compute the tangent vector at snake
  //  * tips.
  //  */
  // unsigned delta_;

  // /*
  //  * Distance threshold for determining overlap with other snakes.
  //  */
  // double overlap_threshold_;

  // /*
  //  * Distance threshold for determining junctions.
  //  */
  // double grouping_distance_threshold_;

  // /*
  //  * Number of points apart to compute the snake branch directions for
  //  * grouping.
  //  */
  // unsigned grouping_delta_;

  // /*
  //  * Angle threshold (in radians) for determining if two snake
  //  * branches are smooth enough to be linked together during grouping
  //  * process.
  //  */
  // double direction_threshold_;

  // /*
  //  * Flag of damping of stretching along z direction. If it is true,
  //  * the stretching forces are reduced if the tangent directions of
  //  * tips are along z direction.
  //  */
  // bool damp_z_;


  DISALLOW_COPY_AND_ASSIGN(Multisnake);
};

} // namespace soax

#endif // SOAX_MULTISNAKE_H_
