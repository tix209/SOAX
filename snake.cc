#include <iomanip>
#include "snake.h"
#include "solver_bank.h"
#include "utility.h"

namespace soax {

SolverBank *Snake::solver_bank_ = NULL;
double Snake::background_ = 0.0;
double Snake::desired_spacing_ = 0.0;
unsigned Snake::max_iterations_ = 0;
double Snake::change_threshold_ = 0.0;
unsigned Snake::check_period_ = 0;
double Snake::minimum_length_ = 0.0;
double Snake::gamma_ = 0.0;
double Snake::external_factor_ = 0.0;
double Snake::stretch_factor_ = 0.0;
int Snake::number_of_sectors_ = 0;
int Snake::radial_near_ = 0;
int Snake::radial_far_ = 0;
unsigned Snake::delta_ = 0;
double Snake::overlap_threshold_ = 0.0;
double Snake::grouping_distance_threshold_ = 0.0;
unsigned Snake::grouping_delta_ = 0;
double Snake::direction_threshold_ = 0.0;
bool Snake::damp_z_ = false;

const double Snake::kBoundary = 0.5;


Snake::Snake(const PointContainer &points, bool is_open, bool is_grouping,
             ImageType::Pointer image,
             VectorImageType::Pointer external_force,
             InterpolatorType::Pointer interpolator,
             VectorInterpolatorType::Pointer vector_interpolator,
             TransformType::Pointer transform) :
    vertices_(points), open_(is_open), grouping_(is_grouping) {
  image_ = image;
  external_force_ = external_force;
  interpolator_ = interpolator;
  vector_interpolator_ = vector_interpolator;
  transform_ = transform;
  viable_ = true;
  initial_state_ = false;
  converged_ = false;
  final_ = false;
  length_ = 0.0;
  spacing_ = 0.0;
  intensity_ = 0.0;
  iterations_ = 0;
  head_tangent_.Fill(0);
  tail_tangent_.Fill(0);
  fixed_head_.Fill(-1.0);
  fixed_tail_.Fill(-1.0);
  head_hooked_snake_ = NULL;
  tail_hooked_snake_ = NULL;
  head_hooked_index_ = 0;
  tail_hooked_index_ = 0;
}


void Snake::Resample() {
  if (!viable_) return;

  if (vertices_.size() < 2) {
    // std::cout << "Snake die: size less than 2!" << std::endl;
    viable_ = false;
    return;
  }

  double spacing = initial_state_ ? 0.25 : desired_spacing_;

  PairContainer sums[kDimension];
  this->UpdateLength(sums);
  unsigned new_size = this->ComputeNewSize(spacing);
  spacing_ = length_ / (new_size - 1);
  this->InterpolateVertices(sums, new_size);

  if (final_ || !open_) {
    viable_ = length_ > minimum_length_;
  } else if (grouping_) {
    viable_ = length_ > grouping_distance_threshold_;
  } else if (!converged_) {
    viable_ = vertices_.size() >= kMinimumEvolvingSize;
  }
}

void Snake::UpdateLength(PairContainer *sums) {
  double current_length = 0.0;
  for (unsigned k = 0; k < kDimension; ++k) {
    sums[k].push_back(std::make_pair(current_length, vertices_.front()[k]));
  }

  for (unsigned i = 1; i < vertices_.size(); ++i) {
    current_length += vertices_[i].EuclideanDistanceTo(vertices_[i-1]);
    for (unsigned k = 0; k < kDimension; ++k) {
      sums[k].push_back(std::make_pair(current_length, vertices_.at(i)[k]));
    }
  }
  length_ = current_length;
}

double Snake::ComputeNewSize(double spacing) const {
  double size = length_ / spacing;
  double diff = length_ - std::floor(size) * spacing;
  if (diff > spacing/2)
    return static_cast<unsigned>(std::ceil(size)) + 1;
  else
    return static_cast<unsigned>(std::floor(size)) + 1;
}


void Snake::InterpolateVertices(const PairContainer *sums,
                                unsigned new_size) {
  PointContainer new_vertices(new_size);

  new_vertices.front() = vertices_.front();
  new_vertices.back() = vertices_.back();

  for (unsigned i = 1; i < new_size-1; ++i) {
    for (unsigned k = 0; k < kDimension; ++k) {
      PairContainer::const_iterator it1, it2;
      it1 = lower_bound(sums[k].begin(), sums[k].end(),
                        std::make_pair(spacing_*i, 0.0));
      it2 = it1--;
      new_vertices.at(i)[k] = it1->second + (it2->second - it1->second) *
          (spacing_*i - it1->first) / (it2->first - it1->first);
    }
  }
  vertices_ = new_vertices;
}


void Snake::Evolve(const SnakeContainer &converged_snakes,
                   unsigned max_iter) {
  unsigned iter = 0;
  // interpolator_->SetInputImage(image_);
  // vector_interpolator_->SetInputImage(external_force_);
  while (iter <= max_iter) {
    if (iterations_ >= max_iterations_)  {
      // std::cout << this << " reaches maximum iterations." << std::endl;
      converged_ = true;
      break;
    }

    if (!(iterations_ % check_period_)) {
      if (this->IsConverged()) break;
      if (iterations_ && initial_state_)
        initial_state_ = false;
    }

    this->CheckSelfIntersection();
    if (!viable_)  break;
    this->HandleHeadOverlap(converged_snakes);
    if (!viable_)  break;
    this->HandleTailOverlap(converged_snakes);
    if (!viable_)  break;
    this->IterateOnce();
    // this->PrintSelf();
    this->Resample();
    iter++;
    if (!viable_)  break;
  }
  this->CheckBodyOverlap(converged_snakes);
}

bool Snake::IsConverged() {
  if (vertices_.size() != last_vertices_.size()) {
    last_vertices_ = vertices_;
    return false;
  }
  for (unsigned i = 0; i < vertices_.size(); ++i) {
    double d = vertices_.at(i).EuclideanDistanceTo(last_vertices_.at(i));

    if (d > change_threshold_) {
      //std::cout << "d: " << d << std::endl;
      last_vertices_ = vertices_;
      return false;
    }
  }
  converged_ = true;
  return true;
}

void Snake::CheckSelfIntersection() {
  if (!open_) return;
  const int min_loop_size = 20;
  if (vertices_.size() <= min_loop_size) return;

  for (PointIterator it1 = vertices_.begin();
       it1 != vertices_.end() - min_loop_size; ++it1) {
    for (PointIterator it2 = it1 + min_loop_size;
         it2 != vertices_.end(); ++it2) {
      double d = (*it1).EuclideanDistanceTo(*it2);
      if (d < 1.0 && !this->TipsStopAtSameLocation()) {
        this->TryInitializeFromPart(vertices_.begin(), it1, true);
        this->TryInitializeFromPart(it2, vertices_.end(), true);
        this->TryInitializeFromPart(it1, it2, false);
        // std::cout << "\nSelf-intersection detected!" << std::endl;
        viable_ = false;
        return;
      }
    }
  }
}

bool Snake::TipsStopAtSameLocation() {
  if (head_hooked_snake_ && tail_hooked_snake_) {
    if (fixed_head_.EuclideanDistanceTo(fixed_tail_) < 1.0)
      return true;
  }
  return false;
}

void Snake::TryInitializeFromPart(PointIterator it1, PointIterator it2,
                                  bool is_open) {
  PointContainer points(it1, it2);
  Snake *s = new Snake(points, is_open, false, image_,
                       external_force_, interpolator_,
                       vector_interpolator_, transform_);
  s->Resample();
  if (s->viable()) {
    subsnakes_.push_back(s);
    // return s;
  } else {
    delete s;
    //return NULL;
  }
}

void Snake::HandleHeadOverlap(const SnakeContainer &converged_snakes) {
  PointIterator first_detach = vertices_.begin();
  PointIterator start = vertices_.begin();

  if (this->HeadIsFixed())
    start += static_cast<int>(overlap_threshold_/spacing_);
  first_detach = this->CheckHeadOverlap(start, converged_snakes);

  if (first_detach != start) {
    if (first_detach == vertices_.end()) {
      // std::cout << "head: total overlap!" << std::endl;
      viable_ = false;
    } else {
      PointIterator last_touch = first_detach - 1;
      this->FindHookedSnakeAndIndex(*last_touch, converged_snakes,
                                    head_hooked_snake_,
                                    head_hooked_index_);
      fixed_head_ = head_hooked_snake_->GetPoint(head_hooked_index_);
      vertices_.erase(vertices_.begin(), first_detach);
      vertices_.push_front(fixed_head_);
      if (!open_)
        open_ = true;
      this->Resample();
    }
  }
}

void Snake::HandleTailOverlap(const SnakeContainer &converged_snakes) {
  PointIterator first_detach = vertices_.end()-1;
  PointIterator start = vertices_.end()-1;

  if (this->TailIsFixed())
    start -= static_cast<int>(overlap_threshold_/spacing_);
  first_detach = this->CheckTailOverlap(start, converged_snakes);


  if (first_detach != start) {
    if (first_detach == vertices_.begin()) {
      // std::cout << "tail: total overlap!" << std::endl;
      viable_ = false;
    } else {
      PointIterator last_touch = first_detach + 1;
      this->FindHookedSnakeAndIndex(*last_touch, converged_snakes,
                                    tail_hooked_snake_,
                                    tail_hooked_index_);
      fixed_tail_ = tail_hooked_snake_->GetPoint(tail_hooked_index_);
      vertices_.erase(last_touch, vertices_.end());
      vertices_.push_back(fixed_tail_);
      if (!open_)
        open_ = true;
      this->Resample();
    }
  }
}

PointIterator Snake::CheckHeadOverlap(
    PointIterator const &start, const SnakeContainer &converged_snakes) {
  PointIterator it = start;
  while (it != vertices_.end() && VertexOverlap(*it, converged_snakes)) {
    ++it;
  }
  return it;
}

PointIterator Snake::CheckTailOverlap(
    PointIterator const &start, const SnakeContainer &converged_snakes) {
  PointIterator it = start;
  while (it != vertices_.begin() && VertexOverlap(*it, converged_snakes)) {
    --it;
  }
  return it;
}

bool Snake::VertexOverlap(const PointType &p,
                          const SnakeContainer &converged_snakes) {
  if (converged_snakes.empty()) return false;
  for (SnakeConstIterator it = converged_snakes.begin();
       it != converged_snakes.end(); ++it) {
    if ((*it)->PassThrough(p, overlap_threshold_))
      return true;
  }
  return false;
}

bool Snake::PassThrough(const PointType &p, double threshold) const {
  for (PointConstIterator it = vertices_.begin();
       it != vertices_.end(); ++it) {
    double dist = p.EuclideanDistanceTo(*it);
    if (dist < threshold)
      return true;
  }
  return false;
}

void Snake::FindHookedSnakeAndIndex(const PointType &p,
                                    const SnakeContainer &converged_snakes,
                                    Snake * &s, unsigned &index) {
  double min_d = kPlusInfinity;
  for (SnakeConstIterator it = converged_snakes.begin();
       it != converged_snakes.end(); ++it) {
    unsigned ind;
    double d = (*it)->FindClosestIndexTo(p, ind);
    if (d < min_d) {
      min_d = d;
      index = ind;
      s = *it;
    }
  }
}


double Snake::FindClosestIndexTo(const PointType &p, unsigned &ind) {
  double min_d = p.EuclideanDistanceTo(vertices_[0]);
  ind = 0;

  for (unsigned i = 1; i < vertices_.size(); ++i) {
    double d = p.EuclideanDistanceTo(vertices_[i]);
    if (d < min_d) {
      min_d = d;
      ind = i;
    }
  }
  return min_d;
}

void Snake::IterateOnce() {
  VectorContainer rhs;
  this->ComputeRHSVector(rhs);
  // this->PrintVectorContainer(rhs);
  for (unsigned direction = 0; direction < kDimension; ++direction) {
    solver_bank_->SolveSystem(rhs, direction, open_);
    this->KeepWithinBounds(direction);
  }

  if (this->HeadIsFixed())
    vertices_.front() = fixed_head_;
  if (this->TailIsFixed())
    vertices_.back() = fixed_tail_;

  iterations_++;
}

void Snake::KeepWithinBounds(unsigned direction) {
  double image_size = image_->GetLargestPossibleRegion().GetSize()[direction];

  for (unsigned i = 0; i < vertices_.size(); ++i) {
    double value = solver_bank_->GetSolution(vertices_.size(), i, open_);
    if (value < kBoundary) {
      vertices_.at(i)[direction] = kBoundary;
    } else if (value > image_size - kBoundary - 1.0) {
      vertices_.at(i)[direction] = image_size - kBoundary - 1.0;
    } else {
      vertices_.at(i)[direction] = value;
    }
  }
}

void Snake::ComputeRHSVector(VectorContainer &rhs) {
  this->AddExternalForce(rhs);
  if (open_)
    this->AddStretchingForce(rhs);
  this->AddVerticesInfo(rhs);
}

void Snake::AddExternalForce(VectorContainer &rhs) {
  for (unsigned i = 0; i < vertices_.size(); ++i) {
    VectorType force = vector_interpolator_->Evaluate(vertices_.at(i));
    force *= external_factor_;
    rhs.push_back(force);
  }
}


void Snake::AddStretchingForce(VectorContainer &rhs) {
  if (!this->HeadIsFixed()) {
    this->UpdateHeadTangent();
    double z_damp = 1;
    if (damp_z_)
      z_damp = exp(-fabs(head_tangent_[2]));

    double head_multiplier = this->ComputeLocalStretch(true);
    rhs.front() += stretch_factor_ * z_damp * head_multiplier * head_tangent_;
  }

  if (!this->TailIsFixed()) {
    this->UpdateTailTangent();
    double z_damp = 1;
    if (damp_z_)
      z_damp = exp(-fabs(tail_tangent_[2]));

    double tail_multiplier = this->ComputeLocalStretch(false);
    rhs.back() += stretch_factor_ * z_damp * tail_multiplier * tail_tangent_;
  }
}

void Snake::AddVerticesInfo(VectorContainer &rhs) {
  for (unsigned i = 0; i < vertices_.size(); ++i) {
    VectorType vector = vertices_.at(i).GetVectorFromOrigin();
    rhs.at(i) += vector * gamma_;
  }
}

void Snake::UpdateHeadTangent() {
  if (vertices_.size() - 1 < delta_)
    head_tangent_ = vertices_.front() - vertices_.back();
  else
    head_tangent_ = vertices_.front() - vertices_.at(delta_);
  head_tangent_.Normalize();
}

void Snake::UpdateTailTangent() {
  if (vertices_.size() - 1 < delta_)
    tail_tangent_ = vertices_.back() - vertices_.front();
  else
    tail_tangent_ = vertices_.back() -
        vertices_.at(vertices_.size() - 1 - delta_);
  tail_tangent_.Normalize();
}


double Snake::ComputeLocalStretch(bool is_head) {
  // PointType &vertex = is_head ? vertices_.front() : vertices_.back();
  // double fg = this->ComputeVertexIntensity(vertex);

  // Good for noisy images such as OCT vessels and microtubules,
  // and actin rings.
  double fg = this->ComputeCircularMeanIntensity(is_head, true);

  if (fg < background_ + kEpsilon)
    return 0.0;

  double bg = this->ComputeCircularMeanIntensity(is_head, false);
  // std::cout << fg << "\t" << bg << std::endl;
  return stretch_factor_ * (1 - bg / fg);
}

double Snake::ComputeCircularMeanIntensity(bool is_head, bool is_fg) {
  const VectorType &normal = is_head ? head_tangent_ : tail_tangent_;
  PointType &vertex = is_head ? vertices_.front() : vertices_.back();
  VectorType radial;
  this->GetStartingRadialDirection(radial, normal, vertex);

  int rnear = is_fg ? 0 : radial_near_;
  int rfar = is_fg ? radial_near_ : radial_far_;

  DataContainer intensities;
  for (int s = 0; s < number_of_sectors_; s++) {
    for (int d = rnear; d < rfar; d++) {
      PointType sample_point;
      this->ComputeSamplePoint(sample_point, vertex, radial, normal, d, s);
      if (this->IsInsideImage(sample_point)) {
        // double intensity = this->ComputeVertexIntensity(sample_point);
        intensities.push_back(interpolator_->Evaluate(sample_point));
      }
    }
  }
  if (intensities.empty()) {
    // std::cout << "empty intensities!" << std::endl;
    // this->PrintSelf();
    return 0; // could be a bug here!!!
  } else {
    return Mean(intensities);
  }
}

void Snake::GetStartingRadialDirection(VectorType &direction,
                                       const VectorType &normal,
                                       const PointType &vertex) {
  unsigned principal_index = this->GetPrincipalIndex(normal);

  PointType other;
  other[(principal_index + 1) % 3] = 0;
  other[(principal_index + 2) % 3] = 0;
  other[principal_index] = (normal[0]*vertex[0] + normal[1]*vertex[1] +
                            normal[2]*vertex[2]) / normal[principal_index];

  direction = other - vertex;
  direction.Normalize();
}


unsigned Snake::GetPrincipalIndex(const VectorType &vec) {
  unsigned ind = 0;
  for (unsigned i = 1; i < kDimension; ++i) {
    if (fabs(vec[i]) > fabs(vec[ind]))
      ind = i;
  }
  return ind;
}

void Snake::ComputeSamplePoint(PointType &point, const PointType &origin,
                               const VectorType &radial,
                               const VectorType &normal,
                               int d, int s) {
  if (d == 0) {
    point = origin;
    return;
  }
  point = origin + static_cast<double>(d) * radial;
  if (s) {
    transform_->SetRotation(normal, 2*kPi*s / number_of_sectors_);
    transform_->SetCenter(origin);
    point = transform_->TransformPoint(point);
  }
}

bool Snake::IsInsideImage(const PointType &point) {
  ImageType::SizeType size = image_->GetLargestPossibleRegion().GetSize();
  for (unsigned i = 0; i < kDimension; ++i) {
    if (point[i] < kBoundary || point[i] >= size[i] - kBoundary)
      return false;
  }
  return true;
}


/*
 * Implementation Notes: CheckBodyOverlap:
 * ---------------------------------------
 * Only try to detect the first overlap part in the middle (not starting from
 * tips) and rely on subsequent execution of this method to detect more
 * overlap part in the body. The reason for detecting first body overlap only
 * is the non-overlap parts will form new snakes which will evolve again.
 */
void Snake::CheckBodyOverlap(const SnakeContainer &converged_snakes) {
  if (!converged_) return;

  bool last_is_overlap = true;
  PointIterator overlap_start = vertices_.end();
  PointIterator overlap_end = vertices_.end();
  // PointType overlap_point;
  for (PointIterator it = vertices_.begin(); it != vertices_.end(); ++it) {
    bool overlap = Snake::VertexOverlap(*it, converged_snakes);
    if (overlap && !last_is_overlap) {
      overlap_start = it;
      //break;
    } else if (!overlap && last_is_overlap && it > overlap_start) {
      overlap_end = it;
      // overlap_point = *it;
      break;
    }
    last_is_overlap = overlap;
  }

  if (overlap_end != vertices_.end()) {
    this->TryInitializeFromPart(overlap_end-1, vertices_.end(), true);
    this->TryInitializeFromPart(vertices_.begin(), overlap_start+1, true);
    viable_ = false;
    // std::cout << "\nbody overlap detected." << std::endl;
  }
}

void Snake::PrintSelf() const {
  std::cout << "\n======== Snake Info ========" << std::endl;
  std::cout << "snake: " << this << std::endl;
  std::cout << "viable: " << viable_ << std::endl;
  std::cout << "open: " << open_ << std::endl;
  std::cout << "converged: " << converged_ << std::endl;
  std::cout << "grouping: " << grouping_ << std::endl;
  std::cout << "iteration: " << iterations_ << std::endl;
  std::cout << "length: " << length_ << std::endl;
  std::cout << "size: " << this->GetSize() << std::endl;
  std::cout << "spacing: " << spacing_ << std::endl;
  std::cout << "fixed head: " << fixed_head_ << std::endl;
  std::cout << "fixed tail: " << fixed_tail_ << std::endl;

  unsigned snake_index(0);
  const unsigned column_width = 15;
  std::cout << "#" << open_ << std::endl;
  for (unsigned j = 0; j != vertices_.size(); ++j) {
    std::cout << snake_index << "\t" << j << "\t";
    std::cout << std::setw(column_width) << this->GetX(j)
              << std::setw(column_width) << this->GetY(j)
              << std::setw(column_width) << this->GetZ(j)
              << std::endl;
  }
}

void Snake::PrintVectorContainer(const VectorContainer &vc) {
  if (vc.empty()) return;
  std::cout << "============ Vector Container  ============" << std::endl;
  for (VectorContainer::const_iterator it = vc.begin();
       it != vc.end(); ++it) {
    std::cout << *it << "\n";
  }
  std::cout << "===========================================" << std::endl;
}

void Snake::UpdateHookedIndices() {
  if (head_hooked_snake_) {
    head_hooked_snake_->AddJunctionIndex(head_hooked_index_);
  }
  if (tail_hooked_snake_) {
    tail_hooked_snake_->AddJunctionIndex(tail_hooked_index_);
  }
}

void Snake::AddJunctionIndex(unsigned index) {
  junction_indices_.insert(index);
}

void Snake::CopySubSnakes(SnakeContainer &c) {
  PointIterator s, e;
  s = e = vertices_.begin();

  if (!junction_indices_.empty()) {
    // add the head sub snake which is not subject to
    // the grouping_distance_threshold
    std::vector<unsigned> indices(junction_indices_.begin(),
                                  junction_indices_.end());
    // for (unsigned i = 0; i < indices.size(); ++i) {
    //   std::cout << indices[i] << std::endl;
    // }
    // std::cout << "end of it" << std::endl;

    s = vertices_.begin() + indices[0];
    PointContainer points(vertices_.begin(), s);
    points.push_back(*s);
    Snake *snake = new Snake(points, true, false, image_,
                             external_force_, interpolator_,
                             vector_interpolator_, transform_);
    snake->Resample();
    if (snake->viable())
      c.push_back(snake);

    // for (IndexSet::iterator it = junction_indices_.begin() + 1;
    //      it != junction_indices_.end(); ++it) {
    //   e = vertices_.begin() + (*it);
    //   PointContainer points(s, e);
    //   // add the junction point too
    //   points.push_back(*e);

    //   Snake *snake = new Snake(points, true, true);
    //   snake->Resample();
    //   if (snake->viable())
    //     c.push_back(snake);
    //   s = e;
    // }
    for (std::vector<unsigned>::iterator it = indices.begin() + 1;
         it != indices.end(); ++it) {
      e = vertices_.begin() + (*it);
      PointContainer points(s, e);
      // add the junction point too
      points.push_back(*e);
      Snake *snake = new Snake(points, true, true, image_,
                               external_force_, interpolator_,
                               vector_interpolator_, transform_);
      snake->Resample();
      if (snake->viable())
        c.push_back(snake);
      s = e;
    }
  }
  // add the tail sub snake which is not subject to
  // the grouping_distance_threshold
  PointContainer points(s, vertices_.end());
  Snake *snake = new Snake(points, true, false, image_,
                           external_force_, interpolator_,
                           vector_interpolator_, transform_);
  snake->Resample();
  if (snake->viable())
    c.push_back(snake);
}

void Snake::EvolveWithTipFixed(unsigned max_iter) {
  unsigned iter = 0;
  fixed_head_ = vertices_.front();
  fixed_tail_ = vertices_.back();

  while (viable_ && iter < max_iter) {
    if (!(iterations_ % check_period_)) {
      if (this->IsConverged())
        break;
    }

    this->IterateOnce();
    this->Resample();
    iter++;
  }
  final_ = true;
  this->Resample();
}

const PointType &Snake::GetTip(bool is_head) const {
  if (is_head)
    return this->GetHead();
  else
    return this->GetTail();
}

} // namespace soax
