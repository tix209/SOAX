#ifndef SOAX_MULTISNAKE_H_
#define SOAX_MULTISNAKE_H_

#include <QObject>
#include "global.h"
#include "snake.h"
#include "junctions.h"


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


  Multisnake();
  ~Multisnake();
  void Reset();
  void ResetContainers();

  /*
   * Load the image and set image_filename_.
   */
  void LoadImage(const std::string &filename);

  /** Load a single file of image sequence.
   */
  void LoadImageSequence(const std::string &filename, int nslices);

  std::string GetImageName(bool suffix = true) const;

  PointType GetImageCenter() const;

  /*
   * Resample and save as an isotropic 16-bit image.
   */
  void SaveAsIsotropicImage(const std::string &filename, double z_spacing);
  void SaveAsIsotropicSequence(const std::string &filename, double z_spacing);

  ImageType::Pointer image() const {return image_;}
  const std::vector<ImageType::Pointer> &image_sequence() const {
    return image_sequence_;
  }
  VectorImageType::Pointer external_force() const {return external_force_;}

  void LoadParameters(const std::string &filename);
  void UpdateSnakeParameters();
  void SaveParameters(const std::string &filename) const;
  void WriteParameters(std::ostream &os) const;

  double intensity_scaling() const  {return intensity_scaling_;}

  /** Set the INTENSITY_SCALING_ to SCALE. SCALE = 0 means
   * automatically scale the maximum intensity to 1.0.
   */
  void set_intensity_scaling(double scale);

  double sigma() const {return sigma_;}
  void set_sigma(double sigma) {sigma_ = sigma;}

  double ridge_threshold() const {return ridge_threshold_;}
  void set_ridge_threshold(double threshold) {
    ridge_threshold_ = threshold;
  }

  unsigned short foreground() const {return foreground_;}
  void set_foreground(unsigned short foreground) {
    foreground_ = foreground;
  }

  unsigned short background() const {return background_;}
  void set_background(unsigned short background) {
    background_ = background;
  }

  bool initialize_z() const {return initialize_z_;}
  void set_initialize_z(bool init_z) {initialize_z_ = init_z;}

  bool is_2d() const {return is_2d_;}

  SolverBank *solver_bank() const {return solver_bank_;}

  /*
   * Invert image intensity by subtracting from the maximum intensity.
   */
  void InvertImageIntensity();

  int GetNumberOfFrames() const {return image_sequence_.size();}

  /*
   * Compute image gradient field for both snake initialization and
   * evolution. If reset is true, the external_force_ is recomputed.
   */
  void ComputeImageGradient(bool reset = true);
  void ComputeImageGradientForSequence(int index);

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
  const std::vector<SnakeContainer> &converged_snakes_sequence() const {
    return converged_snakes_sequence_;
  }

  const std::vector<PointContainer> &junctions_sequence() const {
    return junctions_sequence_;
  }

  void SaveConvergedSnakesAsJFilamentFormat(
      const std::string &filename) const {
    this->SaveJFilamentSnakes(converged_snakes_, filename);
  }

  // void DeformSnakes(QProgressBar * progress_bar = NULL);
  void DeformSnakes();
  void DeformSnakesForSequence();
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

  void LoadSnakesSequence(const std::string &filename);

  void PrintSnakes(const SnakeContainer &snakes) const;

  void SaveSnakes(const SnakeContainer &snakes,
                  const std::string &filename) const;
  void SaveSnakesSequence(const std::string &filename) const;
  // void EvaluateByVertexErrorHausdorffDistance(
  //     const std::string &snake_path, const std::string &filename) const;
  // void EvaluateByFFunction(double threshold, double penalizer,
  //                          int radial_near, int radial_far,
  //                          const std::string &snake_path,
  //                          const std::string &filename) const;

  void PrintGroundTruthLocalSNRValues(int radial_near, int radial_far) const;

  void ComputeRadialOrientation(const PointType &center,
                                double pixel_size,
                                std::ostream & os) const;

  void ComputePointDensityAndIntensity(const PointType &center,
                                       unsigned max_r, double pixel_size,
                                       unsigned type, std::ostream & os) const;

  void ComputeCurvature(int coarse_graining, std::ostream &os) const;
  void ComputeSphericalOrientation(std::ostream &os) const;

  void DeleteSnakes(const SnakeSet &snakes);

  Snake * PopLastInitialSnake();

  void AddInitialSnake(Snake *s) {initial_snakes_.push_back(s);}
  void AddConvergedSnake(Snake *s) {converged_snakes_.push_back(s);}
  void AddSubsnakesToInitialSnakes(Snake *s);

  /*
   * Estimate image SNR using Otsu's method.
   */
  // double ComputeImageSNR(const std::string &binary_filename = "") const;
  // double ComputeImageSNR2(const std::string &filename = "") const;
  // double ComputeForegroundSNR() const;

  // double ComputeGroundTruthFValue(double snr_threshold,
  //                                 double penalizer,
  //                                 int radial_near,
  //                                 int radial_far) const {
  //   return this->ComputeFValue(comparing_snakes1_,
  //                              snr_threshold, penalizer,
  //                              radial_near, radial_far);
  // }

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

  // double ComputeResultSnakesFValue(double snr_threshold,
  //                                  double penalizer,
  //                                  int radial_near,
  //                                  int radial_far) const {
  //   return this->ComputeFValue(converged_snakes_,
  //                              snr_threshold, penalizer,
  //                              radial_near, radial_far);
  // }

  void ComputeResultSnakesVertexErrorHausdorffDistance(
      double &vertex_error, double &hausdorff) const;

  void GenerateSyntheticImage(unsigned foreground,
                              unsigned background,
                              double sigma,
                              const std::string &filename) const;
  void GenerateSyntheticRealImage(double ratio, double sigma,
                                  const std::string &filename) const;

  void GenerateSyntheticTamara(double foreground, double background,
                               const char *filename) const;
  void GenerateSyntheticImageShotNoise(unsigned foreground, unsigned background, unsigned offset,
                                       double scaling, const std::string &filename) const;
  void LoadCurves(const char *filename, double *offset);

  double ComputeDropletMeanIntensity(PointType center, double radius);


 signals:
  void ExtractionProgressed(int value);
  void ExtractionCompleteForFrame(int frame_index);

 private:
  typedef itk::Vector<bool, kDimension> BoolVectorType;
  typedef itk::Image<BoolVectorType, kDimension> BoolVectorImageType;

  void SetImage(int index);

  ImageType::Pointer InterpolateImage(ImageType::Pointer img,
                                      double z_spacing);

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

  void AssignParameters(const std::string &name,
                        const std::string &value);


  void CutSnakes(SnakeContainer &seg);
  void ClearSnakeContainer(SnakeContainer &snakes);
  void ClearSnakeContainerSequence(std::vector<SnakeContainer> &snakes_sequence);


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

  /*
   * Compute F-value for a set of snakes.
   */
  // double ComputeFValue(const SnakeContainer &snakes,
  //                      double threshold, double penalizer,
  //                      int radial_near,int radial_far) const;

  /*
   * Compute a binary image of input image using Otsu's method and
   * return the applied threshold.
   */
  // int ComputeBinaryImage(ImageType::Pointer &img) const;

  /*
   * Compute the intensity statistics of input image using the binary
   * image produced by Otsu's method (see above method)
   */
  // void ComputeForegroundBackgroundStatistics(int threshold,
  //                                            double &fg_mean,
  //                                            double &bg_mean,
  //                                            double &bg_std) const;

  ImageType::PixelType GetMaxImageIntensity() const;

  void ComputeLocalSNRs(const SnakeContainer &snakes,
                        int radial_near, int radial_far,
                        DataContainer &snrs) const;

  void ComputeHistogram(std::vector<unsigned> &hist, int threshold) const;


  double ComputeOtsuThreshold(const std::vector<unsigned> &hist) const;

  void PrintCandidatePoints(BoolVectorImageType::Pointer image,
                            std::ostream &os, unsigned direction) const;




  std::string image_filename_;
  ImageType::Pointer image_;
  std::vector<ImageType::Pointer> image_sequence_;
  VectorImageType::Pointer external_force_;

  InterpolatorType::Pointer interpolator_;
  VectorInterpolatorType::Pointer vector_interpolator_;
  TransformType::Pointer transform_;
  SolverBank *solver_bank_;

  SnakeContainer initial_snakes_;
  SnakeContainer converged_snakes_;
  std::vector<SnakeContainer> converged_snakes_sequence_;
  std::vector<PointContainer> junctions_sequence_;
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
  unsigned short foreground_;
  unsigned short background_;

  /*
   * True if initialize snakes along z axis direction.
   */
  bool initialize_z_;

  /*
   * True if input image is 2d. If true, SOAX behaviour is adapted to
   * 2D. The output SOACs z coordinates is 0.
   */
  bool is_2d_;

  // /*
  //  * True if the intensity of input image needs to be inverted.
  //  */
  // bool invert_intensity_;

  DISALLOW_COPY_AND_ASSIGN(Multisnake);
};

} // namespace soax

#endif // SOAX_MULTISNAKE_H_
