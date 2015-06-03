/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file defines the multiple SOACs (snakes) class for SOAX.
 */


#ifndef MULTISNAKE_H_
#define MULTISNAKE_H_

#include <string>
#include <QObject>  // NOLINT(build/include_order)
#include "./global.h"
#include "./snake.h"
#include "./junctions.h"


class QProgressBar;

namespace soax {

class SolverBank;

/*
 * Multiple open snake class.
 */
class Multisnake : public QObject {
  Q_OBJECT

 public:
  typedef itk::Image<double, kDimension> FloatImageType;

  explicit Multisnake(QObject *parent = 0);
  ~Multisnake();
  void Reset();
  void ResetContainers();
  /*
   * Load the image and set image_filename_.
   */
  void LoadImage(const std::string &filename);

  std::string GetImageName(bool suffix = true) const;

  PointType GetImageCenter() const;

  /**
   * Returns the length of the image diagonal.
   */
  double GetImageDiagonal() const;

  /*
   * Resample and save as an isotropic 16-bit image.
   */
  void SaveAsIsotropicImage(const std::string &filename, double z_spacing);

  ImageType::Pointer image() const {return image_;}
  VectorImageType::Pointer external_force() const {return external_force_;}

  void LoadParameters(const std::string &filename);
  void UpdateSnakeParameters();
  void SaveParameters(const std::string &filename) const;
  std::ostream& WriteParameters(std::ostream &os) const;

  double intensity_scaling() const  {return intensity_scaling_;}

  /** Set the INTENSITY_SCALING_ to SCALE. SCALE = 0 means
   * automatically scale the maximum intensity to 1.0.
   */
  void set_intensity_scaling(double scale);
  double GetIntensityScaling() const;

  double sigma() const {return sigma_;}
  void set_sigma(double sigma) {sigma_ = sigma;}

  double ridge_threshold() const {return ridge_threshold_;}
  void set_ridge_threshold(double threshold) {
    ridge_threshold_ = threshold;
  }

  unsigned foreground() const {return foreground_;}
  void set_foreground(unsigned foreground) {
    foreground_ = foreground;
  }

  unsigned background() const {return background_;}
  void set_background(unsigned background) {
    background_ = background;
  }

  bool initialize_z() const {return initialize_z_;}
  void set_initialize_z(bool init_z) {initialize_z_ = init_z;}

  bool is_2d() const {return is_2d_;}

  SolverBank *solver_bank() const {return solver_bank_;}

  void InvertImageIntensity();

  /*
   * Compute image gradient field for both snake initialization and
   * evolution. If reset is true, the external_force_ is recomputed.
   */
  void ComputeImageGradient(bool reset = true);

  void InitializeSnakes();

  unsigned GetNumberOfInitialSnakes() const {
    return initial_snakes_.size();
  }
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
  const SnakeContainer &converged_snakes() const {
    return converged_snakes_;
  }
  const SnakeContainer &comparing_snakes1() const {
    return comparing_snakes1_;
  }
  const SnakeContainer &comparing_snakes2() const {
    return comparing_snakes2_;
  }

  void SaveConvergedSnakesAsJFilamentFormat(
      const std::string &filename) const {
    this->SaveJFilamentSnakes(converged_snakes_, filename);
  }

  void DeformSnakes();
  void CutSnakesAtTJunctions();
  void GroupSnakes();

  Junctions & junctions() {return junctions_;}
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

  void PrintGroundTruthLocalSNRValues(int radial_near, int radial_far) const;

  void ComputeSphericalOrientation(const PointType &center, double max_r,
                                   double padding, std::ostream &os) const;

  void ComputeRadialOrientation(const PointType &center,
                                double pixel_size,
                                std::ostream &os) const;

  void ComputePointDensityAndIntensity(const PointType &center,
                                       double max_radius, double pixel_size,
                                       std::ostream &os) const;

  void ComputeCurvature(int coarse_graining, double pixel_size,
                        double padding, std::ostream &os) const;

  void ComputeSnakeLength(double pixel_size, std::ostream &os) const;

  void DeleteSnakes(const SnakeSet &snakes);

  Snake * PopLastInitialSnake();

  void AddInitialSnake(Snake *s) {initial_snakes_.push_back(s);}
  void AddConvergedSnake(Snake *s) {converged_snakes_.push_back(s);}
  void AddSubsnakesToInitialSnakes(Snake *s);

  void ComputeGroundTruthLocalSNRs(int radial_near, int radial_far,
                                   DataContainer &snrs) const;
  void ComputeResultSnakesLocalSNRs(int radial_near, int radial_far,
                                    DataContainer &snrs) const;

  double ComputeGroundTruthFValue(const DataContainer &snrs,
                                  double threshold, double penalizer) const;
  double ComputeResultSnakesFValue(const DataContainer &snrs,
                                   double threshold, double penalizer) const;
  double ComputeFValue(const DataContainer &snrs,
                       double threshold, double penalizer) const;

  void ComputeResultSnakesVertexErrorHausdorffDistance(
      double &vertex_error, double &hausdorff) const;

  void GenerateSyntheticImage(unsigned foreground,
                              unsigned background,
                              double sigma,
                              const std::string &filename) const;


 signals:
  void ExtractionProgressed(int value);

 private:
  typedef itk::Vector<bool, kDimension> BoolVectorType;
  typedef itk::Image<BoolVectorType, kDimension> BoolVectorImageType;

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

  static bool IsDarker(Snake *s1, Snake *s2) {
    return s1->ComputeIntensity() < s2->ComputeIntensity();
  }

  static bool IsBrighter(Snake *s1, Snake *s2) {
    return s1->ComputeIntensity() > s2->ComputeIntensity();
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
  void SaveJFilamentSnakes(const SnakeContainer &snakes,
                           const std::string &filename) const;

  void ComputeErrorFromSnakesToComparingSnakes(DataContainer &errors) const;
  void ComputeErrorFromComparingSnakesToSnakes(DataContainer &errors) const;
  double ComputeShortestDistance(const PointType &p,
                                 const SnakeContainer &snakes) const;


  void ComputeRTheta(const PointType &point1, const PointType &point2,
                     const PointType &center, double &r, double &theta) const;

  void ComputeThetaPhi(VectorType vector, double &theta, double &phi) const;

  ImageType::PixelType GetMaxImageIntensity() const;

  void ComputeLocalSNRs(const SnakeContainer &snakes,
                        int radial_near, int radial_far,
                        DataContainer &snrs) const;

  void PrintCandidatePoints(BoolVectorImageType::Pointer image,
                            std::ostream &os, unsigned direction) const;

  bool IsInsideSphere(const PointType &center,
                      double r, const PointType &p) const;

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
  unsigned foreground_;
  unsigned background_;

  /*
   * True if initialize snakes along z axis direction.
   */
  bool initialize_z_;

  /*
   * True if input image is 2d. If true, SOAX behaviour is adapted to
   * 2D. The output SOACs z coordinates is 0.
   */
  bool is_2d_;

  DISALLOW_COPY_AND_ASSIGN(Multisnake);
};

}  // namespace soax

#endif  // MULTISNAKE_H_
