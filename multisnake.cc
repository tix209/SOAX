#include <fstream>
#include <iomanip>
#include <QProgressBar>
#include "multisnake.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkVectorCastImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNormalVariateGenerator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkInvertIntensityImageFilter.h"
#include "solver_bank.h"
#include "utility.h"

namespace soax {

Multisnake::Multisnake() : image_(NULL), external_force_(NULL),
                           intensity_scaling_(0.0), sigma_(1.0),
                           ridge_threshold_(0.01), foreground_(65535),
                           background_(0), initialize_z_(false),
                           is_2d_(false) {
  interpolator_ = InterpolatorType::New();
  vector_interpolator_ = VectorInterpolatorType::New();
  transform_ = TransformType::New();
  solver_bank_ = new SolverBank;
}

Multisnake::~Multisnake() {
  this->ClearSnakeContainer(initial_snakes_);
  this->ClearSnakeContainer(converged_snakes_);
  this->ClearSnakeContainer(comparing_snakes1_);
  this->ClearSnakeContainer(comparing_snakes2_);
  delete solver_bank_;
}

void Multisnake::Reset() {
  this->ClearSnakeContainer(initial_snakes_);
  this->ClearSnakeContainer(converged_snakes_);
  this->ClearSnakeContainer(comparing_snakes1_);
  this->ClearSnakeContainer(comparing_snakes2_);
  junctions_.Reset();
  image_filename_ = "";
  image_ = NULL;
  external_force_ = NULL;
  solver_bank_->Reset();
}

void Multisnake::LoadImage(const std::string &filename) {
  image_filename_ = filename;

  // // Find out the pixel type of the image in file
  // typedef itk::ImageIOBase::IOComponentType  ScalarPixelType;

  // itk::ImageIOBase::Pointer io = itk::ImageIOFactory::CreateImageIO(
  //     filename.c_str(), itk::ImageIOFactory::ReadMode);

  // if( !io ) {
  //   std::cerr << "NO IMAGEIO WAS FOUND" << std::endl;
  //   return;
  // }

  // // Now that we found the appropriate ImageIO class, ask it to
  // // read the meta data from the image file.
  // io->SetFileName(filename);
  // io->ReadImageInformation();
  // ScalarPixelType pixel_type = io->GetComponentType();
  // std::cout << itk::ImageIOBase::GetComponentTypeAsString(pixel_type) << std::endl;
  // unsigned num_dimensions = io->GetNumberOfDimensions();
  // std::cout << num_dimensions << std::endl;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename);

  try {
    reader->Update();
  } catch(itk::ExceptionObject &e) {
    std::cerr << "Exception caught when reading an image!" << std::endl;
    std::cerr << e << std::endl;
  }

  image_ = reader->GetOutput();

  const ImageType::SizeType &size =
      image_->GetLargestPossibleRegion().GetSize();
  std::cout << "Image size: " << size << std::endl;
  if (size[2] < 2) {
    is_2d_ = true;
  }
  interpolator_->SetInputImage(image_);
}

std::string Multisnake::GetImageName(bool suffix) const {
  unsigned last_slash_pos = image_filename_.find_last_of("/\\");
  if (!suffix) {
    unsigned last_dot_pos = image_filename_.find_last_of(".");
    return image_filename_.substr(last_slash_pos+1,
                                  last_dot_pos-last_slash_pos-1);
  } else {
    return image_filename_.substr(last_slash_pos+1);
  }
}

PointType Multisnake::GetImageCenter() const {
  PointType center;
  center.Fill(0.0);
  if (image_) {
    ImageType::SizeType size = image_->GetLargestPossibleRegion().GetSize();
    for (unsigned i = 0; i < kDimension; ++i)
      center[i] = size[i] / 2.0;
  }
  return center;
}

void Multisnake::SaveAsIsotropicImage(const std::string &filename,
                                      double z_spacing) {
  typedef itk::BSplineInterpolateImageFunction<ImageType, double, double>
      InterpolatorType;
  InterpolatorType::Pointer interp = InterpolatorType::New();
  interp->SetSplineOrder(3);
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();
  resampler->SetInterpolator(interp);
  resampler->SetOutputOrigin(image_->GetOrigin());
  resampler->SetOutputDirection(image_->GetDirection());

  // const ImageType::SpacingType &input_spacing = image_->GetSpacing();
  ImageType::SpacingType input_spacing;
  input_spacing.Fill(1.0);
  input_spacing[2] = z_spacing;
  image_->SetSpacing(input_spacing);

  const ImageType::SizeType &input_size =
      image_->GetLargestPossibleRegion().GetSize();
  ImageType::SpacingType output_spacing;
  output_spacing.Fill(1.0);
  ImageType::SizeType output_size;
  for (unsigned i = 0; i < kDimension; ++i) {
    output_size[i] = static_cast<ImageType::SizeValueType>(
        input_size[i] * input_spacing[i]);
  }

  resampler->SetOutputSpacing(output_spacing);
  resampler->SetSize(output_size);
  resampler->SetDefaultPixelValue(0.0);
  resampler->SetInput(image_);

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(resampler->GetOutput());

  try {
    writer->Update();
  } catch( itk::ExceptionObject & exp ) {
    std::cerr << "Exception caught when write an image!" << std::endl;
    std::cerr << exp << std::endl;
  }
}

void Multisnake::LoadParameters(const std::string &filename) {
  std::ifstream infile(filename.c_str());
  if (!infile.is_open()) {
    std::cerr << "LoadParameters: couldn't open file: "
              << infile << std::endl;
    return;
  }

  std::string line;
  while (std::getline(infile, line)) {
    std::stringstream converter;
    converter << line;
    std::string name, value;
    converter >> name >> value;
    this->AssignParameters(name, value);
  }
}

void Multisnake::AssignParameters(const std::string &name,
                                  const std::string &value) {
  if (name == "intensity-scaling") {
    this->set_intensity_scaling(String2Double(value));
    Snake::set_intensity_scaling(intensity_scaling_);
  } else if (name == "smoothing") {
    sigma_ = String2Double(value);
  } else if (name == "grad-diff") {
    ridge_threshold_ = String2Double(value);
  } else if (name == "foreground") {
    foreground_ = String2Unsigned(value);
    Snake::set_foreground(foreground_);
  } else if (name == "background") {
    background_ = String2Unsigned(value);
    Snake::set_background(background_);
  } else if (name == "spacing") {
    Snake::set_desired_spacing(String2Double(value));
  } else if (name == "init-z") {
    initialize_z_ = value == "true";
  } else if (name == "minimum-size") {
    Snake::set_minimum_length(String2Double(value));
  } else if (name == "max-iterations") {
    Snake::set_max_iterations(String2Unsigned(value));
  } else if (name == "change-threshold") {
    Snake::set_change_threshold(String2Double(value));
  } else if (name == "check-period") {
    Snake::set_check_period(String2Unsigned(value));
  } else if (name == "alpha") {
    solver_bank_->set_alpha(String2Double(value));
  } else if (name == "beta") {
    solver_bank_->set_beta(String2Double(value));
  } else if (name == "gamma") {
    solver_bank_->set_gamma(String2Double(value));
    // Snake::set_gamma(solver_bank_->gamma());
  } else if (name == "weight") {
    Snake::set_external_factor(String2Double(value));
  } else if (name == "stretch") {
    Snake::set_stretch_factor(String2Double(value));
  } else if (name == "nsector") {
    Snake::set_number_of_sectors(String2Unsigned(value));
  } else if (name == "radial-near") {
    Snake::set_radial_near(String2Unsigned(value));
  } else if (name == "radial-far") {
    Snake::set_radial_far(String2Unsigned(value));
  } else if (name == "delta") {
    Snake::set_delta(String2Unsigned(value));
  } else if (name == "overlap-threshold") {
    Snake::set_overlap_threshold(String2Double(value));
  } else if (name == "grouping-distance-threshold") {
    Snake::set_grouping_distance_threshold(String2Double(value));
  } else if (name == "grouping-delta") {
    Snake::set_grouping_delta(String2Unsigned(value));
  } else if (name == "direction-threshold") {
    Snake::set_direction_threshold(String2Double(value));
  } else if (name == "damp-z") {
    Snake::set_damp_z(value == "true");
  }
}


void Multisnake::SaveParameters(const std::string &filename) const {
  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    std::cerr << "Couldn't open file: " << outfile << std::endl;
    return;
  }
  this->WriteParameters(outfile);
  outfile.close();
}

void Multisnake::WriteParameters(std::ostream &os) const {
  os << std::boolalpha;
  os << "intensity-scaling\t" << intensity_scaling_ << std::endl;
  os << "smoothing\t" << sigma_ << std::endl;
  os << "grad-diff\t" << ridge_threshold_ << std::endl;
  os << "foreground\t" << foreground_ << std::endl;
  os << "background\t" << background_ << std::endl;
  os << "init-z\t" << initialize_z_ << std::endl;
  os << "spacing \t" << Snake::desired_spacing() << std::endl;
  os << "minimum-size \t" << Snake::minimum_length() << std::endl;
  os << "max-iterations \t" << Snake::max_iterations() << std::endl;
  os << "change-threshold \t" << Snake::change_threshold() << std::endl;
  os << "check-period \t" << Snake::check_period() << std::endl;
  os << "alpha \t" << solver_bank_->alpha() << std::endl;
  os << "beta \t" << solver_bank_->beta() << std::endl;
  os << "gamma \t" << solver_bank_->gamma() << std::endl;
  os << "weight \t" << Snake::external_factor() << std::endl;
  os << "stretch \t" << Snake::stretch_factor() << std::endl;
  os << "nsector \t" << Snake::number_of_sectors() << std::endl;
  os << "radial-near \t" << Snake::radial_near() << std::endl;
  os << "radial-far \t" << Snake::radial_far() << std::endl;
  os << "delta \t" << Snake::delta() << std::endl;
  os << "overlap-threshold \t" << Snake::overlap_threshold() << std::endl;
  os << "grouping-distance-threshold \t"
     << Snake::grouping_distance_threshold() << std::endl;
  os << "grouping-delta \t" << Snake::grouping_delta() << std::endl;
  os << "direction-threshold \t" << Snake::direction_threshold() << std::endl;
  os << "damp-z \t" << Snake::damp_z() << std::endl;
  os << std::noboolalpha;
}

void Multisnake::InvertImageIntensity() {
  ImageType::PixelType maximum = this->GetMaxImageIntensity();
  std::cout << "Previous Maximum intensity: " << maximum << std::endl;

  typedef itk::InvertIntensityImageFilter<ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image_);
  filter->SetMaximum(maximum);

  // typedef itk::RescaleIntensityImageFilter<ImageType, ImageType>
  //     RescalerType;
  // RescalerType::Pointer rescaler = RescalerType::New();
  // rescaler->SetInput(filter->GetOutput());
  // rescaler->SetOutputMinimum(0);
  // rescaler->SetOutputMaximum(255);
  // rescaler->Update();
  // image_ = rescaler->GetOutput();
  // std::cout << "maximum: " << filter->GetMaximum() << std::endl;
  filter->Update();
  image_ = filter->GetOutput();
  interpolator_->SetInputImage(image_);
  maximum = this->GetMaxImageIntensity();
  std::cout << "Current Maximum intensity: " << maximum << std::endl;
}

void Multisnake::ComputeImageGradient(bool reset) {
  if (!reset && external_force_) return;
  external_force_ = NULL;

  typedef itk::Image<double, kDimension> InternalImageType;
  typedef itk::ShiftScaleImageFilter<ImageType,
                                     InternalImageType> ScalerType;
  ScalerType::Pointer scaler = ScalerType::New();
  scaler->SetInput(image_);
  scaler->SetScale(intensity_scaling_);
  scaler->SetShift(0.0);
  scaler->Update();

  if (is_2d_ || sigma_ < 0.01) { // no smoothing
    typedef itk::GradientImageFilter<InternalImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    // filter->SetInput(image_);
    filter->SetInput(scaler->GetOutput());
    typedef itk::VectorCastImageFilter<FilterType::OutputImageType,
                                       VectorImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(filter->GetOutput());
    try {
      caster->Update();
    } catch(itk::ExceptionObject & e) {
      std::cerr << "Exception caught when computing image gradient!\n"
                << e << std::endl;
    }
    external_force_ = caster->GetOutput();
    external_force_->DisconnectPipeline();
  } else {
    typedef itk::GradientRecursiveGaussianImageFilter<
      InternalImageType, VectorImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetSigma(sigma_);
    // filter->SetInput(image_);
    filter->SetInput(scaler->GetOutput());
    try {
      filter->Update();
    } catch(itk::ExceptionObject & e) {
      std::cerr << "Exception caught when computing image gradient!\n"
                << e << std::endl;
    }
    external_force_ = filter->GetOutput();
    external_force_->DisconnectPipeline();
  }
  vector_interpolator_->SetInputImage(external_force_);
}

void Multisnake::InitializeSnakes() {
  this->ClearSnakeContainer(initial_snakes_);
  BoolVectorImageType::Pointer ridge_image =
      InitializeBoolVectorImage();
  this->ScanGradient(ridge_image);

  unsigned num_directions = 2;
  if (!is_2d_ && initialize_z_) num_directions = 3;
  // unsigned num_directions = 1; // for blood vessels
  BoolVectorImageType::Pointer candidate_image =
      InitializeBoolVectorImage();

  for (unsigned d = 0; d < num_directions; ++d) {
    this->GenerateCandidates(ridge_image, candidate_image, d);
  }

  // std::ofstream ofs("x-ridges.txt");
  // this->PrintCandidatePoints(ridge_image, ofs, 0);

  for (unsigned d = 0; d < num_directions; ++d) {
    this->LinkCandidates(candidate_image, d);
  }
  // std::sort(initial_snakes_.begin(), initial_snakes_.end(), IsShorter);
  std::sort(initial_snakes_.begin(), initial_snakes_.end(), IsDarker);
  // this->PrintSnakes(initial_snakes_);
}

Multisnake::BoolVectorImageType::Pointer
Multisnake::InitializeBoolVectorImage() {
  BoolVectorImageType::Pointer image = BoolVectorImageType::New();
  image->SetRegions(image_->GetLargestPossibleRegion());
  image->Allocate();

  BoolVectorImageType::PixelType initial_flag;
  initial_flag.Fill(false);
  image->FillBuffer(initial_flag);

  return image;
}

void Multisnake::ScanGradient(BoolVectorImageType::Pointer ridge_image) {
  typedef itk::ImageRegionIteratorWithIndex<BoolVectorImageType>
      OutputIteratorType;

  OutputIteratorType iter(ridge_image,
                          ridge_image->GetLargestPossibleRegion());

  const unsigned d = is_2d_? 2 : 3;
  for (unsigned i = 0; i < d; ++i) {
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
      VectorImageType::IndexType index = iter.GetIndex();
      VectorImageType::IndexType current_index = index;

      unsigned cnt = 0;
      double grad_comp = external_force_->GetPixel(index)[i];
      // if (invert_intensity_) grad_comp = -grad_comp;

      if (grad_comp < ridge_threshold_) {
        continue;
      } else {
        while (true) {
          current_index[i]++;
          cnt++;

          if (!image_->GetLargestPossibleRegion().IsInside(current_index))
            break;

          double curr_grad_comp = external_force_->GetPixel(current_index)[i];
          // if (invert_intensity_) curr_grad_comp = - curr_grad_comp;

          if (curr_grad_comp > ridge_threshold_) {
            break;
          } else if (curr_grad_comp < -ridge_threshold_) {
            current_index[i] -= cnt/2;
            ridge_image->GetPixel(current_index)[i] = true;
            break;
          }
        }
      }
    }
  }
}

void Multisnake::GenerateCandidates(
    BoolVectorImageType::Pointer ridge_image,
    BoolVectorImageType::Pointer candidate_image, unsigned direction) {
  typedef itk::ImageRegionIteratorWithIndex<BoolVectorImageType>
      OutputIteratorType;
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType>
      ConstIteratorType;

  OutputIteratorType iter(candidate_image,
                          candidate_image->GetLargestPossibleRegion());
  ConstIteratorType iter2(image_, image_->GetLargestPossibleRegion());

  for (iter.GoToBegin(), iter2.GoToBegin(); !iter.IsAtEnd();
       ++iter, ++iter2) {
    if (iter2.Value() > foreground_ || iter2.Value() < background_)
      continue;
    BoolVectorImageType::IndexType index = iter.GetIndex();
    BoolVectorImageType::PixelType ridge_value = ridge_image->GetPixel(index);

    unsigned d = is_2d_ ? 2 : 3;
    if (is_2d_) {
      iter.Value()[direction] = ridge_value[(direction+1) % d];
    } else {
      iter.Value()[direction] = ridge_value[(direction+1) % d] &&
                                ridge_value[(direction+2) % d];
    }
  }
}

void Multisnake::PrintCandidatePoints(
    BoolVectorImageType::Pointer image, std::ostream &os,
    unsigned direction) const {
  typedef itk::ImageRegionIteratorWithIndex<BoolVectorImageType> IteratorType;
  IteratorType it(image, image->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    if (it.Value()[direction])
      os << it.GetIndex() << std::endl;
  }
}

void Multisnake::LinkCandidates(
    BoolVectorImageType::Pointer candidate_image, unsigned direction) {
  ImageType::SizeType size = image_->GetLargestPossibleRegion().GetSize();

  if (is_2d_) {
    for (unsigned c0 = 0; c0 < size[direction]; ++c0) {
      for (unsigned c1 = 0; c1 < size[(direction+1)%2]; ++c1) {
        BoolVectorImageType::IndexType current_index;
        current_index[direction] = c0;
        current_index[(direction+1)%2] = c1;
        current_index[2] = 0;

        if (candidate_image->GetPixel(current_index)[direction])
          LinkFromIndex(candidate_image, current_index, direction);
      }
    }
  } else {
    for (unsigned c0 = 0; c0 < size[direction]; ++c0) {
      for (unsigned c1 = 0; c1 < size[(direction+1)%3]; ++c1) {
        for (unsigned c2 = 0; c2 < size[(direction+2)%3]; ++c2) {
          BoolVectorImageType::IndexType current_index;
          current_index[direction] = c0;
          current_index[(direction+1)%3] = c1;
          current_index[(direction+2)%3] = c2;

          if (candidate_image->GetPixel(current_index)[direction])
            LinkFromIndex(candidate_image, current_index, direction);
        }
      }
    }
  }
}

void Multisnake::LinkFromIndex(
    BoolVectorImageType::Pointer candidate_image,
    BoolVectorImageType::IndexType& index, unsigned direction) {
  PointContainer candidates;

  while (true) {
    // if (!IsInsideImage(point)) break;
    if (!candidate_image->GetLargestPossibleRegion().IsInside(index))
      break;
    PointType point;
    point[0] = index[0];
    point[1] = index[1];
    point[2] = index[2];
    candidates.push_back(point);

    BoolVectorImageType::PixelType pixel_value;
    pixel_value[0] = pixel_value[1] = pixel_value[2] = false;
    candidate_image->SetPixel(index, pixel_value);
    //candidate_image->GetPixel(index)[direction] = false;

    bool next_found = FindNextCandidate(candidate_image, index, direction);
    if (!next_found) break;
  }

  if (candidates.size() > 1) {
    Snake *snake = new Snake(candidates, true, false, image_,
                             external_force_, interpolator_,
                             vector_interpolator_, transform_);
    snake->set_initial_state(true);
    snake->Resample();

    if (snake->viable())
      initial_snakes_.push_back(snake);
    else
      delete snake;
  }
}


bool Multisnake::FindNextCandidate(
    BoolVectorImageType::Pointer candidate_image,
    BoolVectorImageType::IndexType &index, unsigned direction) {
  BoolVectorImageType::IndexType current_index = index;
  index[direction]++;

  if (!candidate_image->GetLargestPossibleRegion().IsInside(index))
    return false;

  if (candidate_image->GetPixel(index)[direction])
    return true;

  if (is_2d_) {
    int d1 = (direction + 1) % 2;
    for (int c1 = current_index[d1] - 1;
         c1 <= current_index[d1] + 1; ++c1) {
      index[d1] = c1;
      index[2] = 0;
      if (candidate_image->GetLargestPossibleRegion().IsInside(index) &&
          candidate_image->GetPixel(index)[direction])
        return true;
    }
  } else {
    int d1 = (direction + 1) % 3;
    int d2 = (direction + 2) % 3;

    for (int c1 = current_index[d1] - 1;
         c1 <= current_index[d1] + 1; ++c1) {
      for (int c2 = current_index[d2] - 1;
           c2 <= current_index[d2] + 1; ++c2) {
        index[d1] = c1;
        index[d2] = c2;
        if (candidate_image->GetLargestPossibleRegion().IsInside(index) &&
            candidate_image->GetPixel(index)[direction])
          return true;
      }
    }
  }
  return false;
}

void Multisnake::DeformSnakes(QProgressBar * progress_bar) {
  unsigned ncompleted = 0;
  std::cout << "Initial # snakes: " << initial_snakes_.size() << std::endl;
  this->ClearSnakeContainer(converged_snakes_);

  // std::ofstream outfile("initial-intensities.txt");
  while (!initial_snakes_.empty()) {
    Snake *snake = initial_snakes_.back();
    initial_snakes_.pop_back();
    // outfile << snake->ComputeIntensity() << std::endl;
    solver_bank_->Reset(false);
    snake->Evolve(solver_bank_, converged_snakes_, kBigNumber, is_2d_);

    if (snake->viable()) {
      converged_snakes_.push_back(snake);
    } else {
      initial_snakes_.insert(initial_snakes_.end(),
                             snake->subsnakes().begin(),
                             snake->subsnakes().end());
      delete snake;
    }

    ncompleted++;
    if (progress_bar) {
      progress_bar->setValue(ncompleted);
    }
    std::cout << "\rRemaining: " << std::setw(6)
              << initial_snakes_.size() << std::flush;
  }
  std::cout << "\nConverged # snakes: " << converged_snakes_.size()
            << std::endl;
}

void Multisnake::CutSnakesAtTJunctions() {
  SnakeContainer segments;
  this->CutSnakes(segments);
  this->ClearSnakeContainer(converged_snakes_);
  converged_snakes_ = segments;
}

void Multisnake::CutSnakes(SnakeContainer &seg) {
  if (converged_snakes_.empty()) return;
  for (SnakeIterator it = converged_snakes_.begin();
       it != converged_snakes_.end(); ++it) {
    (*it)->UpdateHookedIndices();
  }

  for (SnakeIterator it = converged_snakes_.begin();
       it != converged_snakes_.end(); ++it) {
    (*it)->CopySubSnakes(seg);
  }
}

void Multisnake::ClearSnakeContainer(SnakeContainer &snakes) {
  if (snakes.empty()) return;
  for (SnakeContainer::iterator it = snakes.begin();
       it != snakes.end(); ++it) {
    delete *it;
  }
  snakes.clear();
}

void Multisnake::GroupSnakes() {
  junctions_.Initialize(converged_snakes_);
  junctions_.Union();
  junctions_.Configure();
  //  junctions_.PrintTipSets();
  //  junctions_.PrintTips();
  this->LinkSegments(converged_snakes_);
  SnakeIterator it = converged_snakes_.begin();
  while (it != converged_snakes_.end()) {
    solver_bank_->Reset(false);
    (*it)->EvolveWithTipFixed(solver_bank_, 100, is_2d_);
    if ((*it)->viable()) {
      it++;
    } else {
      delete *it;
      // std::cout << "snake  erased!" << std::endl;
      it = converged_snakes_.erase(it);
    }
  }
  this->UpdateJunctions();
}

void Multisnake::LinkSegments(SnakeContainer &seg) {
  SnakeContainer c;
  std::reverse(seg.begin(), seg.end());
  while (!seg.empty()) {
    Snake *segment = seg.back();
    seg.pop_back();
    // log is for detecting loop in linking snake segments
    SnakeSet log;
    PointContainer points;
    bool is_open = true;
    this->LinkFromSegment(segment, seg, log, points, is_open);
    delete segment;
    Snake *s = new Snake(points, is_open, false, image_,
                         external_force_, interpolator_,
                         vector_interpolator_, transform_);
    s->Resample();
    c.push_back(s);
  }
  converged_snakes_ = c;
}

void Multisnake::LinkFromSegment(Snake *s, SnakeContainer &seg,
                                 SnakeSet &log, PointContainer &pc,
                                 bool &is_open) {
  log.insert(s);
  pc = s->vertices();
  SnakeTip * t = junctions_.FindSnakeTip(s, true);
  this->LinkFromSegmentTip(t->neighbor(), pc, is_open, seg, log, true);
  t = junctions_.FindSnakeTip(s, false);
  this->LinkFromSegmentTip(t->neighbor(), pc, is_open, seg, log, false);
}

void Multisnake::LinkFromSegmentTip(SnakeTip *neighbor, PointContainer &pc,
                                    bool &is_open, SnakeContainer &seg,
                                    SnakeSet &log, bool from_head) {
  if (!neighbor) {
    return;
  } else if (log.find(neighbor->snake()) != log.end()) {
    is_open = false;
    return;
  } else {
    this->AddToPointContainer(pc, neighbor->snake(), neighbor->is_head(),
                              from_head);
    SnakeContainer::iterator it = std::find(seg.begin(), seg.end(),
                                            neighbor->snake());
    if (it != seg.end())
      seg.erase(it);
    //seg.remove(neighbor->snake());
    log.insert(neighbor->snake());
    delete neighbor->snake();
    SnakeTip * t = junctions_.FindSnakeTip(neighbor->snake(),
                                           !neighbor->is_head());
    this->LinkFromSegmentTip(t->neighbor(), pc, is_open, seg, log,
                             from_head);
  }
}

void Multisnake::AddToPointContainer(PointContainer &pc, Snake *s,
                                     bool is_head, bool from_head) {
  PointContainer p = s->vertices();

  if (from_head) {
    if (is_head)
      std::reverse(p.begin(), p.end());
    pc.insert(pc.begin(), p.begin(), p.end());
  } else {
    if (!is_head)
      std::reverse(p.begin(), p.end());
    pc.insert(pc.end(), p.begin(), p.end());
  }
}

void Multisnake::UpdateJunctions() {
  if (converged_snakes_.empty())
    junctions_.ClearJunctionPoints();

  const PointContainer &junction_points = junctions_.junction_points();
  PointContainer new_junction_points;
  for (PointConstIterator it = junction_points.begin();
       it != junction_points.end(); ++it) {
    unsigned num_of_close_snake = GetNumberOfSnakesCloseToPoint(*it);
    if (num_of_close_snake > 1)
      new_junction_points.push_back(*it);
  }
  // int size_diff = junction_points.size() - new_junction_points.size();
  // std::cout << "size diff: " << size_diff << std::endl;
  junctions_.set_junction_points(new_junction_points);
}

unsigned Multisnake::GetNumberOfSnakesCloseToPoint(const PointType &p) {
  unsigned num = 0;
  // TODO: improve this dist_threshold
  // const double dist_threshold = grouping_distance_threshold_/2;
  const double dist_threshold = Snake::grouping_distance_threshold();
  for (SnakeContainer::const_iterator it = converged_snakes_.begin();
       it != converged_snakes_.end(); ++it) {
    if ((*it)->PassThrough(p, dist_threshold))
      num++;
  }
  return num;
}

void Multisnake::LoadSnakes(const std::string &filename,
                            SnakeContainer &snakes) {
  this->ClearSnakeContainer(snakes);
  std::ifstream infile(filename.c_str());
  if (!infile.is_open()) {
    std::cerr << "LoadSnakes: couldn't open file: " << filename << std::endl;
    return;
  }
  std::string line, name, value;
  PointContainer points;
  bool is_open = true;
  PointContainer junction_points;

  while (std::getline(infile, line)) {
    if (isalpha(line[0])) {
      std::stringstream converter;
      converter << line;
      converter >> name >> value;
      this->AssignParameters(name, value);
    } else if (line[0] == '#') {
      if (points.size() > 1) {
        Snake *s = new Snake(points, is_open, false, image_, external_force_,
                             interpolator_, vector_interpolator_, transform_);
        snakes.push_back(s);
        s->Resample();
      }
      is_open = (line[1] != '0');
      std::istringstream stream(line);
      std::string dummy;
      double x0, y0, z0, x1, y1, z1;
      stream >> dummy >> x0 >> y0 >> z0 >> x1 >> y1 >> z1;
      points.clear();
    } else if (line[0] == '[') {
      this->LoadPoint(line, junction_points);
    } else {
      std::istringstream stream(line);
      double snake_index, point_index, x, y, z;
      stream >> snake_index >> point_index >> x >> y >> z;
      PointType  snake_point;
      snake_point[0] = x;
      snake_point[1] = y;
      snake_point[2] = z;
      points.push_back(snake_point);
    }
  }
  infile.close();

  if (points.size() > 1) {
    Snake *s = new Snake(points, is_open, false, image_, external_force_,
                         interpolator_, vector_interpolator_, transform_);
    snakes.push_back(s);
    s->Resample();
  }
  junctions_.set_junction_points(junction_points);
}

void Multisnake::LoadJFilamentSnakes(const std::string &filename,
                                     SnakeContainer &snakes) {
  this->ClearSnakeContainer(snakes);
  std::ifstream infile(filename.c_str());
  if (!infile) {
    std::cerr << "Couldn't open file: " << infile << std::endl;
    return;
  }

  std::string line;
  PointContainer points;
  bool is_open = true;
  bool pound  = false;

  while (getline(infile, line)) {
    if (isalpha(line[0])) {
      continue;
    } else if (line[0] == '#') {
      if (points.size() > 1) {
        Snake *s = new Snake(points, is_open, false, image_, external_force_,
                             interpolator_, vector_interpolator_, transform_);
        snakes.push_back(s);
        s->Resample();
      }
      points.clear();
      pound = true;
    } else if (line[0] == '0') {
      if (pound) {
        pound = false;
        continue;
      } else {
        std::istringstream  stream(line);
        double zero, index, x, y, z;
        stream >> zero >> index >> x >> y >> z;
        //std::cout << x << " " << y << " " << z << std::endl;
        PointType  snake_point;
        snake_point[0] = x;
        snake_point[1] = y;
        snake_point[2] = z;
        points.push_back(snake_point);
      }
    }
  }

  infile.close();

  if (points.size() > 1) {
    Snake *s = new Snake(points, is_open, false, image_, external_force_,
                         interpolator_, vector_interpolator_, transform_);
    snakes.push_back(s);
    s->Resample();
  }
}

void Multisnake::LoadPoint(const std::string &s, PointContainer &c) {
  std::istringstream buffer(s);
  char padding;
  double x, y, z;
  buffer >> padding >> x >> padding >> y >> padding >> z >> padding;

  PointType p;
  p[0] = x;
  p[1] = y;
  p[2] = z;
  c.push_back(p);
}

void Multisnake::SaveSnakes(const SnakeContainer &snakes,
                            const std::string &filename) const {
  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    std::cerr << "SaveSnakes: Couldn't open file: " << outfile << std::endl;
    return;
  }

  const unsigned column_width = 16;
  outfile << "image\t" << image_filename_ << std::endl;
  this->WriteParameters(outfile);

  if (snakes.empty()) {
    std::cout << "No snakes to save!" << std::endl;
    return;
  }

  unsigned snake_index = 0;
  // outfile << std::setprecision(12);

  for (SnakeConstIterator it = snakes.begin(); it != snakes.end(); ++it) {
    outfile << "#" << (*it)->open() << std::endl;
    for (unsigned j = 0; j != (*it)->GetSize(); ++j) {
      double intensity = interpolator_->Evaluate((*it)->GetPoint(j));
      outfile << snake_index << std::setw(12) << j
              << std::setw(column_width) << (*it)->GetX(j)
              << std::setw(column_width) << (*it)->GetY(j)
              << std::setw(column_width) << (*it)->GetZ(j)
              << std::setw(20) << intensity << std::endl;
    }
    snake_index++;
  }

  junctions_.PrintJunctionPoints(filename);
  outfile.close();
}

void Multisnake::SaveJFilamentSnakes(const SnakeContainer &snakes,
                                     const std::string &filename) const {
  if (snakes.empty())  {
    std::cerr << "No snakes to save as JFilament snakes!" << std::endl;
    return;
  }

  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    std::cerr << "Couldn't open file: " << outfile << std::endl;
    return;
  }

  outfile << "gamma\t" << solver_bank_->gamma() << std::endl;
  outfile << "weight\t" << Snake::external_factor() << std::endl;
  outfile << "zresolution\t" << 1.0 << std::endl;
  outfile << "background\t" << background_ << std::endl;
  outfile << "alpha\t" << solver_bank_->alpha() << std::endl;
  outfile << "smoothing\t" << sigma_ << std::endl;
  outfile << "stretch\t" << Snake::stretch_factor() << std::endl;
  outfile << "spacing\t" << Snake::desired_spacing() << std::endl;
  outfile << "beta\t" << solver_bank_->beta() << std::endl;
  outfile << "foreground\t" << foreground_ << std::endl;

  for (SnakeConstIterator it = snakes.begin(); it != snakes.end(); ++it) {
    outfile << "#\n0" << std::endl;

    for (unsigned j = 0; j != (*it)->GetSize(); ++j) {
      outfile << "0\t" << j << "\t";
      outfile << (*it)->GetX(j) << "\t"
              << (*it)->GetY(j) << "\t"
              << (*it)->GetZ(j) << std::endl;
    }
  }
  outfile.close();
}

void Multisnake::PrintSnakes(const SnakeContainer &snakes) const {
  if (snakes.empty()) {
    std::cout << "Snake container is empty!" << std::endl;
    return;
  }
  for (SnakeConstIterator it = snakes.begin(); it != snakes.end(); ++it) {
    (*it)->PrintSelf();
  }
}


// void Multisnake::EvaluateByVertexErrorHausdorffDistance(
//     const std::string &snake_path,
//     const std::string &filename) const {
//   std::ofstream outfile;
//   outfile.open(filename.c_str(), std::ios::out | std::ios::app);
//   if (!outfile.is_open()) {
//     std::cerr << "Couldn't open error stat file: " << filename << std::endl;
//     return;
//   }

//   double vertex_error = 20.0;
//   double hausdorff = 60.0;
//   if (!converged_snakes_.empty()) {
//     DataContainer errors1, errors2;
//     this->ComputeErrorFromSnakesToComparingSnakes(errors1);
//     this->ComputeErrorFromComparingSnakesToSnakes(errors2);
//     vertex_error = (Mean(errors1) + Mean(errors2))/2;

//     double max_error1 = Maximum(errors1);
//     double max_error2 = Maximum(errors2);
//     hausdorff = max_error1 > max_error2 ? max_error1 : max_error2;
//   }

//   const unsigned width = 10;
//   outfile << std::setw(width) << ridge_threshold_
//           << std::setw(width) << Snake::stretch_factor()
//           << std::setw(width) << vertex_error
//           << std::setw(width) << hausdorff
//           << "\t" << GetImageName(snake_path) << std::endl;
//   outfile.close();
// }

void Multisnake::ComputeErrorFromSnakesToComparingSnakes(
    DataContainer &errors) const {
  for (SnakeConstIterator iter = converged_snakes_.begin();
       iter != converged_snakes_.end(); ++iter) {
    for (unsigned i = 0; i < (*iter)->GetSize(); ++i) {
      double dist = this->ComputeShortestDistance((*iter)->GetPoint(i),
                                                  comparing_snakes1_);
      errors.push_back(dist);
    }
  }
}

void Multisnake::ComputeErrorFromComparingSnakesToSnakes(
    DataContainer &errors) const {
  for (SnakeConstIterator iter = comparing_snakes1_.begin();
       iter != comparing_snakes1_.end(); ++iter) {
    for (unsigned i = 0; i < (*iter)->GetSize(); ++i) {
      double dist = this->ComputeShortestDistance((*iter)->GetPoint(i),
                                                  converged_snakes_);
      errors.push_back(dist);
    }
  }
}

double Multisnake::ComputeShortestDistance(
    const PointType &p, const SnakeContainer &snakes) const {
  double min_dist = kPlusInfinity;
  for (SnakeConstIterator iter = snakes.begin();
       iter != snakes.end(); ++iter) {
    for (unsigned i = 0; i < (*iter)->GetSize(); ++i) {
      double dist = p.EuclideanDistanceTo((*iter)->GetPoint(i));
      if (dist < min_dist)
        min_dist = dist;
    }
  }
  return min_dist;
}


// void Multisnake::EvaluateByFFunction(double threshold, double penalizer,
//                                      int radial_near, int radial_far,
//                                      const std::string &snake_path,
//                                      const std::string &filename) const {
//   std::ofstream outfile;
//   outfile.open(filename.c_str(), std::ios::out | std::ios::app);
//   if (!outfile.is_open()) {
//     std::cerr << "Couldn't open f-function file: " << filename << std::endl;
//     return;
//   }

//   double fvalue = this->ComputeFValue(converged_snakes_, threshold,
//                                       penalizer, radial_near, radial_far);
//   // if (!converged_snakes_.empty()) {
//   //   // interpolator_->SetInputImage(image_);
//   //   unsigned total_npoints = 0;
//   //   unsigned npoints_low_snr = 0;
//   //   for (SnakeConstIterator it = converged_snakes_.begin();
//   //        it != converged_snakes_.end(); it++) {
//   //     for (unsigned i = 0; i < (*it)->GetSize(); i++) {
//   //       double local_snr = 0.0;
//   //       bool local_bg_defined = (*it)->ComputeLocalSNR(
//   //           i, radial_near, radial_far, local_snr);
//   //       if (local_bg_defined) {
//   //         total_npoints++;
//   //         if (local_snr < threshold) npoints_low_snr++;
//   //       }
//   //     }
//   //   }
//   //   ImageType::SizeType size = image_->GetLargestPossibleRegion().GetSize();
//   //   double normalizer = std::sqrt(static_cast<double>(
//   //       size[1]*size[1]+size[2]*size[2]+size[0]*size[0]));
//   //   fvalue = (penalizer * npoints_low_snr - total_npoints) *
//   //       Snake::desired_spacing() / normalizer;
//   // }

//   const unsigned width = 15;
//   outfile << std::setw(width) << ridge_threshold_
//           << std::setw(width) << Snake::stretch_factor()
//           << std::setw(width) << fvalue
//           << std::setw(width) << threshold
//           << std::setw(width) << penalizer
//           << std::setw(width) << radial_near
//           << std::setw(width) << radial_far
//           << std::endl;
//   outfile.close();
// }

void Multisnake::PrintGroundTruthLocalSNRValues(int radial_near,
                                                int radial_far) {
  if (comparing_snakes1_.empty()) return;

  DataContainer snrs;
  for (SnakeConstIterator it = comparing_snakes1_.begin();
       it != comparing_snakes1_.end(); ++it) {
    for (unsigned i = 0; i < (*it)->GetSize(); i++) {
      double local_snr = 0.0;
      bool local_bg_defined = (*it)->ComputeLocalSNRAtIndex(
          i, radial_near, radial_far, local_snr);
      if (local_bg_defined) {
        snrs.push_back(local_snr);
      }
    }
  }

  std::cout << "========= Local SNR Info =========" << std::endl;
  std::cout << "Min: " << Minimum(snrs) << "\tMax: " << Maximum(snrs)
            << std::endl;
  double mean_snr = Mean(snrs);
  std::cout << "Mean: " << mean_snr << "\t Median: " << Median(snrs)
            << "\t Std: " << StandardDeviation(snrs, mean_snr) << std::endl;
}


void Multisnake::ComputeRadialOrientation(const PointType &center,
                                          const std::string &filename) const {
  if (converged_snakes_.empty()) return;

  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    std::cerr << "Cannot open file for radial orientation results."
              << std::endl;
    return;
  }

  const unsigned width = 15;
  outfile << "Image file\t" << image_filename_ << "\n"
          << "Image center\t" << center << "\n"
          << std::setw(width) << "radius (um)"
          << std::setw(width) << "theta" << std::endl;

  const double spacing = 0.166; // droplet image

  for (SnakeConstIterator it = converged_snakes_.begin();
       it != converged_snakes_.end(); ++it) {
    for (unsigned i = 0; i < (*it)->GetSize() - 1; ++i) {
      double r, theta;
      this->ComputeRTheta((*it)->GetPoint(i), (*it)->GetPoint(i+1),
                          center, r, theta);
      outfile << std::setw(width) << r * spacing
              << std::setw(width) << theta << std::endl;
    }
  }

  outfile.close();
}

void Multisnake::ComputeRTheta(const PointType &point1,
                               const PointType &point2,
                               const PointType &center,
                               double &r, double &theta) const {
  PointType mid;
  mid.SetToMidPoint(point1, point2);
  r = mid.EuclideanDistanceTo(center);

  VectorType vector = point1 - point2;
  VectorType radial = mid - center;
  vector.Normalize();
  radial.Normalize();
  theta = std::acos(vector * radial) * 180 / kPi;
  theta = theta > 90.0 ? (180.0 - theta) : theta;
}

void Multisnake::ComputePointDensityAndIntensity(const PointType &center,
                                                 unsigned max_r,
                                                 double pixel_size,
                                                 std::ostream & os) const {
  if (converged_snakes_.empty()) return;

  const unsigned width = 16;
  os << "Image file\t" << image_filename_
     << "\nImage center\t" << center
     << "\nMax radius\t" << max_r << "\n"
     << std::setw(width) << "radius (um)"
     << std::setw(width) << "soac density"
     << std::setw(width) << "soac intensity"
     << std::setw(width) << "voxel intensity" << std::endl;

  DataContainer snake_intensities(max_r, 0.0);
  DataContainer voxel_intensities(max_r, 0.0);
  std::vector<unsigned> snaxel_counts(max_r, 0);
  std::vector<unsigned> voxel_counts(max_r, 0);

  // snaxel counts and snake intensities
  for (SnakeConstIterator it = converged_snakes_.begin();
       it != converged_snakes_.end(); ++it) {
    for (unsigned i = 0; i < (*it)->GetSize(); ++i) {
      const PointType &p = (*it)->GetPoint(i);
      unsigned r = static_cast<unsigned>(center.EuclideanDistanceTo(p));
      if (r < max_r) {
        snaxel_counts[r]++;
        snake_intensities[r] += interpolator_->Evaluate(p);
      }
    }
  }

  // voxel counts and voxel intensities
  itk::ImageRegionConstIteratorWithIndex<ImageType>
      it(image_, image_->GetLargestPossibleRegion());
  it.GoToBegin();
  while (!it.IsAtEnd()) {
    ImageType::IndexType index = it.GetIndex();
    PointType p;
    p[0] = index[0];
    p[1] = index[1];
    p[2] = index[2];
    unsigned r = static_cast<unsigned>(center.EuclideanDistanceTo(p));
    if (r < max_r) {
      voxel_counts[r]++;
      voxel_intensities[r] += it.Value();
    }
    ++it;
  }

  // compute density and average intensity
  for (unsigned i = 0; i < max_r; ++i) {
    if (snaxel_counts[i] > 0)
      snake_intensities[i] /= snaxel_counts[i];

    voxel_intensities[i] /= voxel_counts[i];
    double r = (i+1) * pixel_size;
    double density = snaxel_counts[i] / (4 * kPi * r * r);
    os << std::setw(width) << r
       << std::setw(width) << density
       << std::setw(width) << snake_intensities[i]
       << std::setw(width) << voxel_intensities[i] << std::endl;
  }
}

void Multisnake::ComputeCurvature(int coarse_graining,
                                  const std::string &filename) const {
  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    std::cerr << "Cannot open file for snake point density results."
              << std::endl;
    return;
  }

  outfile << "Image file\t" << image_filename_ << "\n"
          << "Coarse graining\t" << coarse_graining << std::endl;

  for (SnakeConstIterator it = converged_snakes_.begin();
       it != converged_snakes_.end(); ++it) {
    const int end_index = (*it)->GetSize() - 2*coarse_graining;
    for (int i = 0; i < end_index; i += coarse_graining) {
      PointType p0 = (*it)->GetPoint(i);
      PointType p1 = (*it)->GetPoint(i + coarse_graining);
      PointType p2 = (*it)->GetPoint(i + 2*coarse_graining);
      VectorType vec1 = p1 - p0;
      vec1.Normalize();
      VectorType vec2 = p2 - p1;
      vec2.Normalize();
      double length = coarse_graining * (*it)->spacing();
      double curvature = (vec1 - vec2).GetNorm() / length;
      outfile << curvature << std::endl;
    }
  }
  outfile.close();
}

void Multisnake::ComputeSphericalOrientation(
    const std::string &filename) const {
  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    std::cerr << "Cannot open file for snake point density results."
              << std::endl;
    return;
  }

  const unsigned width = 15;
  outfile << "Image file\t" << image_filename_ << "\n"
          << std::setw(width) << "theta"
          << std::setw(width) << "phi" << std::endl;
  for (SnakeConstIterator it = converged_snakes_.begin();
       it != converged_snakes_.end(); ++it) {
    for (unsigned i = 0; i < (*it)->GetSize() - 1; ++i) {
      VectorType vector = (*it)->GetPoint(i) - (*it)->GetPoint(i+1);
      double theta, phi;
      this->ComputeThetaPhi(vector, theta, phi);
      outfile << std::setw(width) << theta
              << std::setw(width) << phi << std::endl;
    }
  }
  outfile.close();
}

void Multisnake::ComputeThetaPhi(VectorType vector,
                                 double &theta, double &phi) const {
  // phi is (-pi/2, +pi/2]
  // theta is [0, pi)
  const double r = vector.GetNorm();

  if (std::abs(vector[0]) < kEpsilon && std::abs(vector[1]) < kEpsilon) {
    // x = y = 0
    phi = 0;
    theta = 0;
  } else if (std::abs(vector[0]) < kEpsilon) {
    // x = 0, y != 0
    if (vector[1] < -kEpsilon)
      vector = -vector;

    phi = 90;
    theta = std::acos(vector[2]/r) * 180 / kPi;
  } else {
    // x != 0
    if (vector[0] < -kEpsilon)
      vector = -vector;

    theta = std::acos(vector[2]/r) * 180 / kPi;
    phi = std::atan(vector[1] / vector[0]) * 180 / kPi;
  }
}

void Multisnake::DeleteSnakes(const SnakeSet &snakes) {
  for (SnakeSet::const_iterator it = snakes.begin();
       it != snakes.end(); ++it) {
    SnakeIterator it2 = std::find(converged_snakes_.begin(),
                                  converged_snakes_.end(), *it);
    if (it2 != converged_snakes_.end()) {
      converged_snakes_.erase(it2);
    } else {
      std::cout << "Couldn't locate snake " << *it << " in converged snakes."
                << std::endl;
    }
  }
}

Snake * Multisnake::PopLastInitialSnake() {
  Snake *s = initial_snakes_.back();
  initial_snakes_.pop_back();
  return s;
}


void Multisnake::AddSubsnakesToInitialSnakes(Snake *s) {
  initial_snakes_.insert(initial_snakes_.end(),
                         s->subsnakes().begin(), s->subsnakes().end());
}

double Multisnake::ComputeImageSNR2(const std::string &filename) const {
  typedef itk::OtsuMultipleThresholdsImageFilter<ImageType,
                                                 ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image_);
  filter->SetNumberOfThresholds(2);
  filter->Update();
  if (!filename.empty()) {
    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(filename);
    writer->SetInput(filter->GetOutput());
    writer->Update();
  }

  FilterType::ThresholdVectorType thresholds = filter->GetThresholds();
  std::cout << "thresholds: " << std::endl;
  for (unsigned i = 0; i < thresholds.size(); ++i) {
    std::cout << thresholds[i] << std::endl;
  }

  // Compute intensity statistics of foreground and background voxels
  // foreground voxels have intensity [thresholds[1], max_intensity]
  // background voxels have intensity [thresholds[0], thresholds[1])
  DataContainer fgs, bgs;
  typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
  ConstIteratorType it(image_, image_->GetLargestPossibleRegion());
  it.GoToBegin();
  while (!it.IsAtEnd()) {
    if (it.Get() >= thresholds[1]) { // foreground
      fgs.push_back(it.Get());
    } else if (it.Get() >= thresholds[0]) { // background
      bgs.push_back(it.Get());
    }
    ++it;
  }

  double fg_mean = Mean(fgs);
  double bg_mean = Mean(bgs);
  double bg_std = StandardDeviation(bgs, bg_mean);
  if (bg_std < kEpsilon) {
    std::cerr << "Background intensity std is zero!" << std::endl;
    return 0.0;
  }
  return (fg_mean - bg_mean) / bg_std;
}

double Multisnake::ComputeImageSNR(const std::string &binary_filename) const {
  ImageType::Pointer img = NULL;
  int threshold = this->ComputeBinaryImage(img);
  std::cout << "Otsu threshold: " << threshold << std::endl;

  if (!binary_filename.empty()) {
    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    // std::string::size_type slash_pos = image_filename_.find_last_of("/\\");
    // std::string::size_type dot_pos = image_filename_.find_last_of(".");
    // std::string extracted_name = image_filename_.substr(
    //     slash_pos+1, dot_pos-slash_pos-1);
    // std::cout << "extracted name: " << extracted_name << std::endl;
    // writer->SetFileName(extracted_name + "_binary.mha");
    writer->SetFileName(binary_filename);
    writer->SetInput(img);
    writer->Update();
  }

  double fg_mean(0.0), bg_mean(0.0), bg_std(0.0);
  this->ComputeForegroundBackgroundStatistics(threshold, fg_mean,
                                              bg_mean, bg_std);
  if (bg_std < kEpsilon) {
    std::cerr << "Background intensity std is zero!" << std::endl;
    return 0.0;
  }
  return (fg_mean - bg_mean) / bg_std;
}

int Multisnake::ComputeBinaryImage(ImageType::Pointer &img) const {
  typedef itk::OtsuThresholdImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer otsu_filter = FilterType::New();
  otsu_filter->SetInput(image_);
  otsu_filter->SetOutsideValue(255);
  otsu_filter->SetInsideValue(0);
  otsu_filter->Update();
  img = otsu_filter->GetOutput();
  return otsu_filter->GetThreshold();
}

void Multisnake::ComputeForegroundBackgroundStatistics(
    int threshold, double &fg_mean, double &bg_mean, double &bg_std) const {
  DataContainer fgs, bgs;
  typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
  ConstIteratorType it(image_, image_->GetLargestPossibleRegion());
  it.GoToBegin();
  while (!it.IsAtEnd()) {
    if (it.Get() > threshold) { // foreground
      fgs.push_back(it.Get());
    } else { // background
      bgs.push_back(it.Get());
    }
    ++it;
  }
  fg_mean = Mean(fgs);
  bg_mean = Mean(bgs);
  bg_std = StandardDeviation(bgs, bg_mean);
  // std::cout << "fg mean: " << fg_mean << "\tbg mean: " << bg_mean
  //           << "\tbg std: " << bg_std << std::endl;
}

double Multisnake::ComputeForegroundSNR() const {
  const unsigned max_intensity_levels = 1 << 16;
  std::vector<unsigned> histogram(max_intensity_levels, 0);
  int intensity_threshold = 0;
  // only consider those voxels with intensity greater than
  // "intensity_threshold"
  this->ComputeHistogram(histogram, intensity_threshold);
  int otsu_threshold1 = this->ComputeOtsuThreshold(histogram);
  std::cout << "first Otsu's threshold = " << otsu_threshold1 << std::endl;
  // reset the histogram
  std::fill(histogram.begin(), histogram.end(), 0);
  this->ComputeHistogram(histogram, otsu_threshold1);
  int otsu_threshold2 = this->ComputeOtsuThreshold(histogram);
  std::cout << "Second Otsu's threshold = " << otsu_threshold2 << std::endl;

  DataContainer fgs, bgs;
  typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
  ConstIteratorType it(image_, image_->GetLargestPossibleRegion());
  it.GoToBegin();
  while (!it.IsAtEnd()) {
    if (it.Get() >= otsu_threshold1) {
      if (it.Get() >= otsu_threshold2) { // foreground
        fgs.push_back(it.Get());
      } else { // background
        bgs.push_back(it.Get());
      }
    }
    ++it;
  }
  double fg_mean = Mean(fgs);
  double bg_mean = Mean(bgs);
  double bg_std = StandardDeviation(bgs, bg_mean);
  if (bg_std < kEpsilon) {
    std::cerr << "Background intensity std is zero!" << std::endl;
    return 0.0;
  }
  return (fg_mean - bg_mean) / bg_std;
}

void Multisnake::ComputeHistogram(std::vector<unsigned> &hist,
                                  int threshold) const {
  typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
  ConstIteratorType it(image_, image_->GetLargestPossibleRegion());
  it.GoToBegin();
  while (!it.IsAtEnd()) {
    if (it.Get() >= threshold) {
      hist[it.Get()]++;
    }
    ++it;
  }
}

double Multisnake::ComputeOtsuThreshold(
    const std::vector<unsigned> &hist) const {
  unsigned total_nvoxels = 0;
  double total_intensities = 0.0;
  for (unsigned i = 0; i < hist.size(); ++i) {
    total_nvoxels += hist[i];
    total_intensities += i * hist[i];
  }

  double var_max = 0.0;
  unsigned threshold = 0;

  double sum_bg = 0.0;
  unsigned weight_bg = 0;
  unsigned weight_fg = 0;
  for (unsigned i = 0; i < hist.size(); ++i) {
    weight_bg += hist[i];
    if (!weight_bg) continue;

    weight_fg = total_nvoxels - weight_bg;
    if (!weight_fg) break;

    sum_bg += static_cast<double>(i * hist[i]);
    double mean_bg = sum_bg / weight_bg;
    double mean_fg = (total_intensities - sum_bg) / weight_fg;
    double var = static_cast<double>(
        weight_bg * weight_fg * (mean_bg - mean_fg) * (mean_bg - mean_fg));

    if (var > var_max) {
      var_max = var;
      threshold = i;
    }
  }

  return threshold;
}

// double Multisnake::ComputeFValue(const SnakeContainer &snakes,
//                                  double threshold, double penalizer,
//                                  int radial_near, int radial_far) const {
//   double fvalue = 0.0;
//   if (!snakes.empty()) {
//     unsigned total_npoints = 0;
//     unsigned low_snr_npoints = 0;
//     for (SnakeConstIterator it = snakes.begin(); it != snakes.end(); it++) {
//       for (unsigned i = 0; i < (*it)->GetSize(); i++) {
//         double local_snr = 0.0;

//         bool local_bg_defined = (*it)->ComputeLocalSNR(
//             i, radial_near, radial_far, local_snr);
//         if (local_bg_defined) {
//           total_npoints++;
//           if (local_snr < threshold) low_snr_npoints++;
//         }
//       }
//     }
//     // std::cout << "c: " << penalizer << std::endl;
//     // std::cout << "Low SNR points: " << low_snr_npoints << std::endl;
//     fvalue = (penalizer * low_snr_npoints - total_npoints) *
//         Snake::desired_spacing();
//   }
//   return fvalue;
// }
double Multisnake::ComputeGroundTruthFValue(const DataContainer &snrs,
                                            double threshold,
                                            double penalizer) const {
  return this->ComputeFValue(snrs, threshold, penalizer) /
      comparing_snakes1_.size();
}

double Multisnake::ComputeResultSnakesFValue(const DataContainer &snrs,
                                             double threshold,
                                             double penalizer) const {
  return this->ComputeFValue(snrs, threshold, penalizer) /
      converged_snakes_.size();
}

double Multisnake::ComputeFValue(const DataContainer &snrs,
                                 double threshold, double penalizer) const {
  if (snrs.empty()) return 0.0;
  unsigned low_snr_npoints = 0;
  for (DataContainer::const_iterator it = snrs.begin(); it != snrs.end();
       ++it) {
    if (*it < threshold) low_snr_npoints++;
  }
  return penalizer * low_snr_npoints - snrs.size();
}

void Multisnake::ComputeGroundTruthLocalSNRs(int radial_near, int radial_far,
                                             DataContainer &snrs) const {
  return ComputeLocalSNRs(comparing_snakes1_, radial_near, radial_far, snrs);
}

void Multisnake::ComputeResultSnakesLocalSNRs(int radial_near, int radial_far,
                                              DataContainer &snrs) const {
  return ComputeLocalSNRs(converged_snakes_, radial_near, radial_far, snrs);
}

void Multisnake::ComputeLocalSNRs(const SnakeContainer &snakes,
                                  int radial_near, int radial_far,
                                  DataContainer &snrs) const {
  if (snakes.empty()) return;
  for (SnakeConstIterator it = snakes.begin(); it != snakes.end(); it++) {
    for (unsigned i = 0; i < (*it)->GetSize(); i++) {
      double local_snr = 0.0;
      bool local_bg_defined = (*it)->ComputeLocalSNRAtIndex(
          i, radial_near, radial_far, local_snr);

      if (local_bg_defined)
        snrs.push_back(local_snr);
      // else
      //   std::cerr << "Local background undefined!" << std::endl;
    }
  }
}

void Multisnake::ComputeResultSnakesVertexErrorHausdorffDistance(
    double &vertex_error, double &hausdorff) const {
  if (converged_snakes_.empty() || comparing_snakes1_.empty()) return;

  DataContainer errors1, errors2;

  this->ComputeErrorFromSnakesToComparingSnakes(errors1);
  this->ComputeErrorFromComparingSnakesToSnakes(errors2);
  vertex_error = (Mean(errors1) + Mean(errors2))/2;
  // std::cout << "ve: " << vertex_error << std::endl;
  double max_error1 = Maximum(errors1);
  double max_error2 = Maximum(errors2);
  hausdorff = max_error1 > max_error2 ? max_error1 : max_error2;
  // std::cout << "hausdorff: " << hausdorff << std::endl;
}

void Multisnake::GenerateSyntheticImage(unsigned foreground,
                                        unsigned background,
                                        double sigma,
                                        const std::string &filename) const {
  if (!image_ || comparing_snakes1_.empty()) return;
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(image_->GetLargestPossibleRegion());
  image->Allocate();
  image->FillBuffer(0);

  // Assign centerline intensities
  if (foreground) {
    for (SnakeConstIterator it = comparing_snakes1_.begin();
         it != comparing_snakes1_.end(); ++it) {
      for (unsigned i = 0; i < (*it)->GetSize(); ++i) {
        ImageType::IndexType index;
        image_->TransformPhysicalPointToIndex((*it)->GetPoint(i), index);
        if (!image->GetPixel(index)) { // only set once
          image->SetPixel(index, foreground);
        }
      }
    }

    // apply PSF
    typedef itk::RecursiveGaussianImageFilter<ImageType,
                                              ImageType>  FilterType;
    FilterType::Pointer filter_x = FilterType::New();
    FilterType::Pointer filter_y = FilterType::New();
    FilterType::Pointer filter_z = FilterType::New();
    filter_x->SetDirection(0);
    filter_y->SetDirection(1);
    filter_z->SetDirection(2);
    filter_x->SetOrder(FilterType::ZeroOrder);
    filter_y->SetOrder(FilterType::ZeroOrder);
    filter_z->SetOrder(FilterType::ZeroOrder);
    filter_x->SetInput(image);
    filter_y->SetInput(filter_x->GetOutput());
    filter_z->SetInput(filter_y->GetOutput());
    const double psf = 2.0;
    filter_x->SetSigma(psf);
    filter_y->SetSigma(psf);
    filter_z->SetSigma(psf);

    // rescale the intensity back to foreground
    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType>
        RescalerType;
    RescalerType::Pointer rescaler = RescalerType::New();
    rescaler->SetInput(filter_z->GetOutput());
    rescaler->SetOutputMinimum(0.0);
    rescaler->SetOutputMaximum(foreground);
    rescaler->Update();
    image = rescaler->GetOutput();
  }

  // add Gaussian noise with u=background and std=sigma
  typedef itk::ImageRegionIterator<ImageType> IteratorType;
  IteratorType iter(image, image->GetLargestPossibleRegion());

  typedef itk::Statistics::NormalVariateGenerator GeneratorType;
  GeneratorType::Pointer generator = GeneratorType::New();
  generator->Initialize(2003);

  for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
    iter.Value() += background + sigma * generator->GetVariate();
  }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(image);
  try {
    writer->Update();
  } catch(itk::ExceptionObject &e) {
    std::cerr << "Snake::SaveImageUShort: Exception caught!\n"
              << e << std::endl;
  }
  return;
}

void Multisnake::set_intensity_scaling(double scale) {
  if (scale > kEpsilon)
    intensity_scaling_ = scale;
  else
    intensity_scaling_ = 1.0 / this->GetMaxImageIntensity();
}

ImageType::PixelType Multisnake::GetMaxImageIntensity() const {
  if (!image_) return 0;
  typedef itk::MinimumMaximumImageCalculator<ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetImage(image_);
  filter->ComputeMaximum();
  return filter->GetMaximum();
}

} // namespace soax
