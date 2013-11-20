#include <fstream>
#include <iomanip>
#include <QProgressBar>
#include "multisnake.h"
#include "itkImageFileReader.h"
#include "itkShiftScaleImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkVectorCastImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "solver_bank.h"
#include "utility.h"

namespace soax {

Multisnake::Multisnake() :
    image_(NULL), external_force_(NULL), intensity_scaling_(0.0),
    intensity_scaled_(false), sigma_(0.0), ridge_threshold_(0.0),
    foreground_(0.0), background_(0.0), initialize_z_(false) {
  interpolator_ = InterpolatorType::New();
  vector_interpolator_ = VectorInterpolatorType::New();
  transform_ = TransformType::New();
  solver_bank_ = new SolverBank;
  Snake::set_solver_bank(solver_bank_);
}

Multisnake::~Multisnake() {
  delete solver_bank_;
}

void Multisnake::LoadImage(const std::string &filename) {
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename);
  image_ = reader->GetOutput();
  try {
    reader->Update();
  } catch(itk::ExceptionObject &e) {
    std::cerr << "Exception caught when reading an image!" << std::endl;
    std::cerr << e << std::endl;
  }
  image_filename_ = filename;
  const ImageType::SizeType &size =
      image_->GetLargestPossibleRegion().GetSize();
  std::cout << "Image size: " << size << std::endl;
  std::cout << image_filename_ << std::endl;
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
  Snake::set_background(background_ * intensity_scaling_);
}

void Multisnake::AssignParameters(const std::string &name,
                                  const std::string &value) {
  if (name == "intensity-scaling") {
    intensity_scaling_ = String2Double(value);
  } else if (name == "smoothing") {
    sigma_ = String2Double(value);
  } else if (name == "grad-diff") {
    ridge_threshold_ = String2Double(value);
  } else if (name == "foreground") {
    foreground_ = String2Double(value);
  } else if (name == "background") {
    background_ = String2Double(value);
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
    Snake::set_gamma(solver_bank_->gamma());
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


// void Multisnake::UpdateSnakeParameters() {
//   Snake::set_background(background_ * intensity_scaling_);
//   Snake::set_desired_spacing(desired_spacing_);
//   Snake::set_minimum_length(minimum_length_);
//   Snake::set_max_iterations(max_iterations_);
//   Snake::set_change_threshold(change_threshold_);
//   Snake::set_check_period(check_period_);
//   Snake::set_gamma(gamma_);
//   Snake::set_external_factor(external_factor_);
//   Snake::set_stretch_factor(stretch_factor_);
//   Snake::set_number_of_sectors(number_of_sectors_);
//   Snake::set_radial_near(radial_near_);
//   Snake::set_radial_far(radial_far_);
//   Snake::set_delta(delta_);
//   Snake::set_overlap_threshold(overlap_threshold_);
//   Snake::set_grouping_distance_threshold(grouping_distance_threshold_);
//   Snake::set_grouping_delta(grouping_delta_);
//   Snake::set_direction_threshold(direction_threshold_);
//   Snake::set_damp_z(damp_z_);
//   Snake::set_solver_bank(solver_bank_);
//   solver_bank_->set_alpha(alpha_);
//   solver_bank_->set_beta(beta_);
//   solver_bank_->set_gamma(gamma_);
// }

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

// void Multisnake::PrintParameters() const {
//   std::cout << "============ Current Parameters ============" << std::endl;
//   std::cout << std::boolalpha;
//   std::cout << "intensity-scaling: " << intensity_scaling_ << std::endl;
//   std::cout << "smoothing: " << sigma_ << std::endl;
//   std::cout << "grad-diff: " << ridge_threshold_ << std::endl;
//   std::cout << "foreground: " << foreground_ << std::endl;
//   std::cout << "background: " << background_ << std::endl;
//   std::cout << "init-z: " << initialize_z_ << std::endl;
//   std::cout << "spacing: " << desired_spacing_ << std::endl;
//   std::cout << "minimum-size: " << minimum_length_ << std::endl;
//   std::cout << "max-iterations: " << max_iterations_ << std::endl;
//   std::cout << "change-threshold: " << change_threshold_ << std::endl;
//   std::cout << "check-period: " << check_period_ << std::endl;
//   std::cout << "alpha: " << alpha_ << std::endl;
//   std::cout << "beta: " << beta_ << std::endl;
//   std::cout << "gamma: " << gamma_ << std::endl;
//   std::cout << "weight: " << external_factor_ << std::endl;
//   std::cout << "stretch: " << stretch_factor_ << std::endl;
//   std::cout << "nsector: " << number_of_sectors_ << std::endl;
//   std::cout << "radial-near: " << radial_near_ << std::endl;
//   std::cout << "radial-far: " << radial_far_ << std::endl;
//   std::cout << "delta: " << delta_ << std::endl;
//   std::cout << "overlap-threshold: " << overlap_threshold_ << std::endl;
//   std::cout << "grouping-distance-threshold: " << grouping_distance_threshold_
//             << std::endl;
//   std::cout << "grouping-delta: " << grouping_delta_ << std::endl;
//   std::cout << "direction-threshold: " << direction_threshold_ << std::endl;
//   std::cout << "damp-z: " << damp_z_ << std::endl;
//   std::cout << "============================================" << std::endl;
//   std::cout << std::noboolalpha;
// }

void Multisnake::ScaleImageIntensity() {
  typedef itk::ShiftScaleImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image_);
  filter->SetScale(intensity_scaling_);
  filter->SetShift(0.0);
  image_ = filter->GetOutput();
  try {
    filter->Update();
  } catch (itk::ExceptionObject &e) {
    std::cerr << "Exception caught when scaling image intensity!\n"
              << e << std::endl;
  }
  intensity_scaled_ = true;
  interpolator_->SetInputImage(image_);
}

void Multisnake::ComputeImageGradient() {
  if (sigma_ < 0.01) { // no smoothing
    typedef itk::GradientImageFilter<ImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(image_);
    typedef itk::VectorCastImageFilter<FilterType::OutputImageType,
                                       VectorImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(filter->GetOutput());
    try {
      caster->Update();
    } catch( itk::ExceptionObject & e ) {
      std::cerr << "Exception caught when computing image gradient!\n"
                << e << std::endl;
    }
    external_force_ = caster->GetOutput();
    external_force_->DisconnectPipeline();
  } else {
    typedef itk::GradientRecursiveGaussianImageFilter<
      ImageType, VectorImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetSigma(sigma_);
    filter->SetInput(image_);
    try {
      filter->Update();
    } catch( itk::ExceptionObject & e ) {
      std::cerr << "Exception caught when computing image gradient!\n"
                << e << std::endl;
    }
    external_force_ = filter->GetOutput();
    external_force_->DisconnectPipeline();
  }
  vector_interpolator_->SetInputImage(external_force_);
}

void Multisnake::InitializeSnakes() {
  BoolVectorImageType::Pointer ridge_image =
      InitializeBoolVectorImage();
  ScanGradient(ridge_image);

  unsigned num_directions = initialize_z_ ? 3 : 2;
  // unsigned num_directions = 1; // for blood vessels
  BoolVectorImageType::Pointer candidate_image =
      InitializeBoolVectorImage();

  for (unsigned d = 0; d < num_directions; ++d) {
    GenerateCandidates(ridge_image, candidate_image, d);
  }

  for (unsigned d = 0; d < num_directions; ++d) {
    LinkCandidates(candidate_image, d);
  }
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

  for (unsigned i = 0; i < kDimension; ++i) {
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
      VectorImageType::IndexType index = iter.GetIndex();
      VectorImageType::IndexType current_index = index;

      unsigned cnt = 0;
      if (external_force_->GetPixel(index)[i] < ridge_threshold_) {
        continue;
      } else {
        while (true) {
          current_index[i]++;
          cnt++;

          if (!image_->GetLargestPossibleRegion().IsInside(current_index))
            break;

          if (external_force_->GetPixel(current_index)[i] >
              ridge_threshold_) {
            break;
          } else if (external_force_->GetPixel(current_index)[i] <
                     -ridge_threshold_) {
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

  const double foreground = foreground_ * intensity_scaling_;
  const double background = background_ * intensity_scaling_;

  for (iter.GoToBegin(), iter2.GoToBegin(); !iter.IsAtEnd();
       ++iter, ++iter2) {
    if (iter2.Value() > foreground || iter2.Value() < background)
      continue;
    BoolVectorImageType::IndexType index = iter.GetIndex();

    iter.Value()[direction] =
        ridge_image->GetPixel(index)[(direction+1) % kDimension]
        && ridge_image->GetPixel(index)[(direction+2) % kDimension];
  }
}

void Multisnake::LinkCandidates(
    BoolVectorImageType::Pointer candidate_image, unsigned direction) {
  ImageType::SizeType size = image_->GetLargestPossibleRegion().GetSize();

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
  }
}


bool Multisnake::FindNextCandidate(
    BoolVectorImageType::Pointer candidate_image,
    BoolVectorImageType::IndexType &index, unsigned direction) {
  BoolVectorImageType::IndexType current_index = index;
  index[direction]++;
  int d1 = (direction + 1) % 3;
  int d2 = (direction + 2) % 3;

  if (candidate_image->GetPixel(index)[direction])
    return true;

  for (int c1 = current_index[d1] - 1;
       c1 <= current_index[d1] + 1; ++c1) {
    for (int c2 = current_index[d2] - 1;
         c2 <= current_index[d2] + 1; ++c2) {
      index[d1] = c1;
      index[d2] = c2;
      if (candidate_image->GetPixel(index)[direction])
        return true;
    }
  }
  return false;
}

void Multisnake::SortSnakesOnLength(SnakeContainer &snakes) {
  std::sort(snakes.begin(), snakes.end(), IsShorter);
}


void Multisnake::SaveSnakes(const SnakeContainer &snakes,
                            const std::string &filename) const {
  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    std::cerr << "Couldn't open file: " << outfile << std::endl;
    return;
  }

  const unsigned column_width = 20;
  outfile << "image\t" << image_filename_ << std::endl;
  WriteParameters(outfile);

  if (snakes.empty()) {
    std::cout << "No snakes to save!" << std::endl;
    return;
  }

  unsigned snake_index = 0;
  outfile << std::setprecision(12);

  for (SnakeConstIterator it = snakes.begin(); it != snakes.end(); ++it) {
    outfile << "#" << (*it)->open() << std::endl;
    for (unsigned j = 0; j != (*it)->GetSize(); ++j) {
      double intensity = interpolator_->Evaluate((*it)->GetPoint(j));
      if (intensity_scaled_)
        intensity /= intensity_scaling_;
      outfile << snake_index << "\t" << j << "\t";
      outfile << std::setw(column_width) << (*it)->GetX(j)
              << std::setw(column_width) << (*it)->GetY(j)
              << std::setw(column_width) << (*it)->GetZ(j)
              << std::setw(column_width) << intensity
              << std::endl;
    }
    snake_index++;
  }

  junctions_.PrintJunctionPoints(filename);
  outfile.close();
}

void Multisnake::DeformSnakes(QProgressBar * progress_bar) {
  unsigned ncompleted = 0;
  std::cout << "Initial # snakes: " << initial_snakes_.size() << std::endl;

  while (!initial_snakes_.empty()) {
    Snake *snake = initial_snakes_.back();
    initial_snakes_.pop_back();
    // if (initial_snakes_.size() == 100)
    // snake->PrintSelf();

    snake->Evolve(converged_snakes_, kBigNumber);

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
    (*it)->EvolveWithTipFixed(100);
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

  PointContainer &junction_points = junctions_.junction_points();
  PointContainer new_junction_points;
  for (PointContainer::iterator it = junction_points.begin();
       it != junction_points.end(); ++it) {
    unsigned num_of_close_snake = GetNumberOfSnakesCloseToPoint(*it);
    if (num_of_close_snake > 1)
      new_junction_points.push_back(*it);
  }
  // int size_diff = junction_points.size() - new_junction_points.size();
  // std::cout << "size diff: " << size_diff << std::endl;
  junction_points = new_junction_points;
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
  snakes.clear();
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

void Multisnake::EvaluateByVertexErrorHausdorffDistance(
    const std::string &snake_path,
    const std::string &filename) const {
  std::ofstream outfile;
  outfile.open(filename.c_str(), std::ios::out | std::ios::app);
  if (!outfile.is_open()) {
    std::cerr << "Couldn't open error stat file: " << filename << std::endl;
    return;
  }

  double vertex_error = static_cast<double>(
      image_->GetLargestPossibleRegion().GetSize()[0]);
  double hausdorff = static_cast<double>(
      image_->GetLargestPossibleRegion().GetSize()[0]);
  if (!converged_snakes_.empty()) {
    DataContainer errors1, errors2;
    this->ComputeErrorFromSnakesToComparingSnakes(errors1);
    this->ComputeErrorFromComparingSnakesToSnakes(errors2);
    vertex_error = (Mean(errors1) + Mean(errors2))/2;

    double max_error1 = Maximum(errors1);
    double max_error2 = Maximum(errors2);
    hausdorff = max_error1 > max_error2 ? max_error1 : max_error2;
  }

  const unsigned width = 10;
  outfile << std::setw(width) << ridge_threshold_
          << std::setw(width) << Snake::stretch_factor()
          << std::setw(width) << vertex_error
          << std::setw(width) << hausdorff
          << "\t" << GetImageName(snake_path) << std::endl;
  outfile.close();
}

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


void Multisnake::EvaluateByFFunction(double threshold, double penalizer,
                                     const std::string &snake_path,
                                     const std::string &filename) const {
}

} // namespace soax
