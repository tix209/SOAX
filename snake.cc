/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file implements the SOAC (snake) class for SOAX.
 */


#include <iomanip>
#include "./snake.h"
#include "./solver_bank.h"
#include "./utility.h"

namespace soax {

double Snake::intensity_scaling_ = 0.0;
unsigned Snake::foreground_ = 65535;
unsigned Snake::background_ = 0;
double Snake::desired_spacing_ = 1.0;
double Snake::minimum_length_ = 10.0;
unsigned Snake::max_iterations_ = 10000;
double Snake::change_threshold_ = 0.1;
unsigned Snake::check_period_ = 100;
unsigned Snake::iterations_per_press_ = 100;
double Snake::external_factor_ = 1.0;
double Snake::stretch_factor_ = 0.2;
int Snake::number_of_sectors_ = 8;
int Snake::radial_near_ = 4;
int Snake::radial_far_ = 8;
int Snake::radial_save_foreground_ = 0;
unsigned Snake::delta_ = 4;
double Snake::overlap_threshold_ = 1.0;
double Snake::grouping_distance_threshold_ = 4.0;
unsigned Snake::grouping_delta_ = 8.0;
double Snake::direction_threshold_ = 2.1;
bool Snake::damp_z_ = false;
double Snake::z_spacing_ = 2.88;
const double Snake::kBoundary = 0.5;


Snake::Snake(const PointContainer &points, bool is_open, bool is_grouping,
             ImageType::Pointer image,
             VectorImageType::Pointer external_force,
             InterpolatorType::Pointer interpolator,
             VectorInterpolatorType::Pointer vector_interpolator,
             TransformType::Pointer transform) :
    open_(is_open), grouping_(is_grouping) {
  vertices_ = points;
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
  if (length_ < spacing) {
    viable_ = false;
    return;
  }
  unsigned new_size = this->ComputeNewSize(spacing);
  spacing_ = length_ / (new_size - 1);
  this->InterpolateVertices(sums, new_size);

  if (final_) {
    viable_ = length_ > minimum_length_;
  } else if (grouping_) {
    viable_ = length_ > grouping_distance_threshold_;
  } else {
    viable_ = vertices_.size() >= kMinimumEvolvingSize;
  }

  // if (final_ || !open_) {
  //   viable_ = length_ > minimum_length_;
  // } else if (grouping_) {
  //   viable_ = length_ > grouping_distance_threshold_;
  // } else if (!converged_) {
  //   viable_ = vertices_.size() >= kMinimumEvolvingSize;
  // }
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


void Snake::Evolve(SolverBank *solver, const SnakeContainer &converged_snakes,
                   unsigned max_iter, unsigned dim, const std::vector<std::vector<IndexPairContainer > >  &converged_snakes_grid) {
  unsigned iter = 0;

  
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
    this->HandleHeadOverlap(converged_snakes, converged_snakes_grid);
    if (!viable_)  break;
    this->HandleTailOverlap(converged_snakes, converged_snakes_grid);
    if (!viable_)  break;
    this->IterateOnce(solver, dim);
    this->Resample();
    iter++;
    if (!viable_)  break;    
  }
  this->CheckBodyOverlap(converged_snakes, converged_snakes_grid);
}

bool Snake::IsConverged() {
  if (vertices_.size() != last_vertices_.size()) {
    last_vertices_ = vertices_;
    return false;
  }
  for (unsigned i = 0; i < vertices_.size(); ++i) {
    double d = vertices_.at(i).EuclideanDistanceTo(last_vertices_.at(i));

    if (d > change_threshold_) {
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
  } else {
    delete s;
  }
}

void Snake::HandleHeadOverlap(const SnakeContainer &converged_snakes, const std::vector<std::vector<IndexPairContainer > > &converged_snakes_grid) {
  PointIterator first_detach = vertices_.begin();
  PointIterator start = vertices_.begin();

  if (this->HeadIsFixed())
    start += static_cast<int>(overlap_threshold_/spacing_);

  //first vertex in current snake that does not have an overlap
  first_detach = this->CheckHeadOverlap(start, converged_snakes, converged_snakes_grid);
  
  
  
  // if there is an overlap  
  if (first_detach != start) {
    if (first_detach == vertices_.end()) {
      // std::cout << "head: total overlap!" << std::endl;
      viable_ = false;
    } else {
      
      // finds closest point to last_touch (last vertex that has an overlap) on a converged_snake
      PointIterator last_touch = first_detach - 1;
      this->FindHookedSnakeAndIndex(*last_touch, converged_snakes,
                                    head_hooked_snake_,
                                    head_hooked_index_, converged_snakes_grid);
                                    
      // point on converged snake closest to last_touch
      fixed_head_ = head_hooked_snake_->GetPoint(head_hooked_index_);
      
      
      // connects last_touch to point on converged snake
      vertices_.erase(vertices_.begin(), first_detach);
      vertices_.push_front(fixed_head_);
      if (!open_)
        open_ = true;
      this->Resample();
    }
  }
}

void Snake::HandleTailOverlap(const SnakeContainer &converged_snakes, const std::vector<std::vector<IndexPairContainer > > &converged_snakes_grid) {
  PointIterator first_detach = vertices_.end()-1;
  PointIterator start = vertices_.end()-1;

  if (this->TailIsFixed())
    start -= static_cast<int>(overlap_threshold_/spacing_);
  first_detach = this->CheckTailOverlap(start, converged_snakes, converged_snakes_grid);


  if (first_detach != start) {
    if (first_detach == vertices_.begin()) {
      // std::cout << "tail: total overlap!" << std::endl;
      viable_ = false;
    } else {
        
      // finds closest point to last_touch on a converged_snake
      PointIterator last_touch = first_detach + 1;
      
      this->FindHookedSnakeAndIndex(*last_touch, converged_snakes,
                                    tail_hooked_snake_,
                                    tail_hooked_index_, converged_snakes_grid);
                                    
      // point on converged snake closest to last_touch
      fixed_tail_ = tail_hooked_snake_->GetPoint(tail_hooked_index_);
      
      // connects last_touch to point on converged snake
      vertices_.erase(last_touch, vertices_.end());
      vertices_.push_back(fixed_tail_);
      if (!open_)
        open_ = true;
      this->Resample();
    }
  }
}

// starting from head vertex checks if each vertex in the called snake overlaps with another converged snake
// returns the first vertex that has no overlap
PointIterator Snake::CheckHeadOverlap(
    PointIterator const &start, const SnakeContainer &converged_snakes, const std::vector<std::vector<IndexPairContainer > >  &converged_snakes_grid) {
        
  int loopCount = 0;
  PointIterator it = start;
  while (it != vertices_.end() && VertexOverlap(*it, converged_snakes, converged_snakes_grid)) {
    loopCount++;
    ++it;
  }
  return it;
}

PointIterator Snake::CheckTailOverlap(
    PointIterator const &start, const SnakeContainer &converged_snakes, const std::vector<std::vector<IndexPairContainer > > &converged_snakes_grid) {
  PointIterator it = start;
  while (it != vertices_.begin() && VertexOverlap(*it, converged_snakes, converged_snakes_grid)) {
    --it;
  }
  return it;
}

// checks if there is an overlap between a point and any of the converged snakes within a threshold distance
bool Snake::VertexOverlap(const PointType &p,
                          const SnakeContainer &converged_snakes, const std::vector<std::vector<IndexPairContainer > > &converged_snakes_grid) {
  if (converged_snakes.empty()) return false;
  /*for (SnakeConstIterator it = converged_snakes.begin();
       it != converged_snakes.end(); ++it) {
    if ((*it)->PassThrough(p, overlap_threshold_))
      return true;
  }*/
  
  int org_x_grid = (int)(p[0] / overlap_threshold_);
  int org_y_grid = (int)(p[1] / overlap_threshold_);

  // bounds check for whether point p is in grid, puts particles outside grid into most extreme existing grid (should be okay, just means that more particles are being compared)
  if(org_x_grid < 0)
  {
      org_x_grid = 0;
  }
  else if(org_x_grid >= converged_snakes_grid.size())
  {
      org_x_grid = converged_snakes_grid.size() - 1;
  }
  
  if(org_y_grid < 0)
  {
      org_y_grid = 0;
  }
  else if(org_y_grid >= converged_snakes_grid[0].size())
  {
      org_y_grid = converged_snakes_grid[0].size() - 1;
  }
  
  IndexPairContainer indexes_of_points_in_grid = converged_snakes_grid[org_x_grid][org_y_grid];

  bool test = false;
  //#pragma omp parallel for
  for(unsigned i = 0; i < indexes_of_points_in_grid.size(); ++i)
  {
      if(test == false)
      {
          PointType point_in_grid = converged_snakes[std::get<0>(indexes_of_points_in_grid[i])]->GetPoint(std::get<1>(indexes_of_points_in_grid[i]));
          double dist = p.EuclideanDistanceTo(point_in_grid);
          if (dist < overlap_threshold_)
          {
              //#pragma omp critical
              {
                test = true;
              }
          }
      }
  }
  
  return test;
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

// loops through the converged snakes to find the closest point on each converged snake to p
// then finds the closest overall point to p on a converged snake
// returns the converged snake and the index of the vertex which is closest to p
void Snake::FindHookedSnakeAndIndex(const PointType &p,
                                    const SnakeContainer &converged_snakes,
                                    Snake * &s, unsigned &index, const std::vector<std::vector<IndexPairContainer > >  &converged_snakes_grid) {
  double min_d = kPlusInfinity;
  
  
  /*for (SnakeConstIterator it = converged_snakes.begin();
       it != converged_snakes.end(); ++it) {
    unsigned ind;
    double d = (*it)->FindClosestIndexTo(p, ind);
    if (d < min_d) {
      min_d = d;
      index = ind;
      s = *it;
    }
  }*/
  
  int org_x_grid = (int)(p[0] / overlap_threshold_);
  int org_y_grid = (int)(p[1] / overlap_threshold_);

  // bounds check for whether point p is in grid, puts particles outside grid into most extreme existing grid (should be okay, just means that more particles are being compared)
  if(org_x_grid < 0)
  {
      org_x_grid = 0;
  }
  else if(org_x_grid >= converged_snakes_grid.size())
  {
      org_x_grid = converged_snakes_grid.size() - 1;
  }
  
  if(org_y_grid < 0)
  {
      org_y_grid = 0;
  }
  else if(org_y_grid >= converged_snakes_grid[0].size())
  {
      org_y_grid = converged_snakes_grid[0].size() - 1;
  }
  
  //int org_z_grid = (int)(p[2] / overlap_threshold_);
  
  IndexPairContainer indexes_of_points_in_grid = converged_snakes_grid[org_x_grid][org_y_grid];
  
  //https://stackoverflow.com/questions/28258590/using-openmp-to-get-the-index-of-minimum-element-parallelly
  //#pragma omp parallel
  {
      int index_local = 0;
      double min_local = min_d;
      Snake *s_local;
      
      //#pragma omp for nowait
      for (unsigned i = 0; i < indexes_of_points_in_grid.size(); ++i) {
        //unsigned ind;
        //double d = converged_snakes[i]->FindClosestIndexTo(p, ind);
        
        int s_local_tmp = std::get<0>(indexes_of_points_in_grid[i]);
        int index_local_tmp = std::get<1>(indexes_of_points_in_grid[i]);
        
        PointType point_in_grid = converged_snakes[s_local_tmp]->GetPoint(index_local_tmp);
        double d = p.EuclideanDistanceTo(point_in_grid);
        
        if (d < min_local) {
          min_local = d;
          index_local = index_local_tmp;
          s_local = converged_snakes[s_local_tmp];
        }
      }
      //#pragma omp critical
      {
          if(min_local < min_d) {
              min_d = min_local;
              index = index_local;
              s = s_local;
          }
      }
  }
      
}

// finds the index of the called snake that is closest to p
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

void Snake::IterateOnce(SolverBank *solver, unsigned dim) {
  VectorContainer rhs;
  this->ComputeRHSVector(solver->gamma(), rhs, dim);
  
  //std::cout << "dim: " << dim << ", verticesSize: " << vertices_.size() << std::endl;
  
  for (unsigned d = 0; d < dim; ++d) {
    solver->SolveSystem(rhs, d, open_);
    for (unsigned i = 0; i < vertices_.size(); ++i) {
      double value = solver->GetSolution(vertices_.size(), i, open_);
      vertices_.at(i)[d] = value;
    }
  }

  if (this->HeadIsFixed())
    vertices_.front() = fixed_head_;
  if (this->TailIsFixed())
    vertices_.back() = fixed_tail_;

  iterations_++;
}

void Snake::ComputeRHSVector(double gamma, VectorContainer &rhs, unsigned dim) {
  this->AddVerticesInfo(gamma, rhs);
  this->AddExternalForce(rhs, dim);
  if (open_)
    this->AddStretchingForce(rhs, dim);
}

void Snake::AddExternalForce(VectorContainer &rhs, unsigned dim) {
  for (unsigned i = 0; i < vertices_.size(); ++i) {
    if (IsInsideImage(vertices_.at(i), dim)) {
      {
      rhs.at(i) += external_factor_ *
                   vector_interpolator_->Evaluate(vertices_.at(i));
      }
    }
  }
}


void Snake::AddStretchingForce(VectorContainer &rhs, unsigned dim) {
  this->UpdateHeadTangent();
  if (IsInsideImage(vertices_.front(), dim) && !this->HeadIsFixed()) {
    double z_damp = 1;
    if (damp_z_)
      z_damp = exp(-fabs(head_tangent_[2]));

    double head_multiplier = this->ComputeLocalStretch(0, dim);
    rhs.front() += stretch_factor_ * z_damp * head_multiplier * head_tangent_;
  }

  this->UpdateTailTangent();
  if (IsInsideImage(vertices_.back(), dim) && !this->TailIsFixed()) {
    double z_damp = 1;
    if (damp_z_)
      z_damp = exp(-fabs(tail_tangent_[2]));

    double tail_multiplier = this->ComputeLocalStretch(
        vertices_.size() - 1, dim);
    rhs.back() += stretch_factor_ * z_damp * tail_multiplier * tail_tangent_;
  }
}

void Snake::AddVerticesInfo(double gamma, VectorContainer &rhs) {
  for (unsigned i = 0; i < vertices_.size(); ++i) {
    VectorType vector = vertices_.at(i).GetVectorFromOrigin();
    rhs.push_back(vector * gamma);
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


/* Note that scaling the intensity is not necessary here because they
 * are cancelled out. */
double Snake::ComputeLocalStretch(unsigned index, unsigned dim) {
  double fg = InterpolateImageIntensity(vertices_[index]);
  if (fg < background_ + kEpsilon || fg > foreground_)
    return 0.0;

  double bg = 0.0;
  if (dim == 2)
    bg = this->ComputeBackgroundMeanIntensity2d(index);
  else
    bg = this->ComputeBackgroundMeanIntensity(index);

  if (fg < bg || bg < 0.0)
    return 0.0;
  return 1.0 - bg / fg;
}

double Snake::ComputeBackgroundMeanIntensity(unsigned index) const {
  const PointType &vertex = vertices_[index];
  const VectorType &normal = this->ComputeUnitTangentVector(index);

  assert(std::fabs(normal.GetNorm() - 1.0) < 1e-9);
  VectorType z_axis;
  z_axis[0] = 0.0;
  z_axis[1] = 0.0;
  z_axis[2] = 1.0;
  VectorType projection = z_axis - normal[2] * normal;
  VectorType long_axis;
  long_axis[0] = 1.0;
  long_axis[1] = 0.0;
  long_axis[2] = 0.0;
  VectorType short_axis;
  short_axis[0] = 0.0;
  short_axis[1] = 1.0;
  short_axis[2] = 0.0;
  double projection_length = projection.GetNorm();
  if (projection_length > kEpsilon) {
    long_axis = projection / projection_length;
    short_axis = itk::CrossProduct(long_axis, normal);
    short_axis.Normalize();
  }
  DataContainer bgs;
  const double angle_step = 2 * kPi / number_of_sectors_;
  for (int r = radial_near_; r < radial_far_; r++) {
    for (int s = 0; s < number_of_sectors_; s++) {
      double angle = s * angle_step;
      VectorType v = static_cast<double>(r) * (std::cos(angle) * long_axis +
                                               std::sin(angle) * short_axis);
      v[2] *= z_spacing_;
      PointType p = vertex + v;
      if (IsInsideImage(p)) {
        double intensity = InterpolateImageIntensity(p);
        if (intensity > background_)
          bgs.push_back(intensity);
      }
    }
  }

  if (bgs.empty()) {
    return -1.0;  // return a negative value intentionally
  } else {
    return Mean(bgs);
  }
}

double Snake::ComputeBackgroundMeanIntensity2d(unsigned index) const {
  const VectorType &normal = this->ComputeUnitTangentVector(index);
  PointType vertex = vertices_[index];
  DataContainer bgs;

  for (int d = radial_near_; d < radial_far_; d++) {
    PointType pod;
    pod[0] = this->ComputePodX(vertex[0], normal, d, true);
    pod[1] = this->ComputePodY(vertex[1], normal, d, false);
    pod[2] = vertex[2];

    if (IsInsideImage(pod, 2)) {
      double intensity = InterpolateImageIntensity(pod);
      bgs.push_back(intensity);
    }

    pod[0] = this->ComputePodX(vertex[0], normal, d, false);
    pod[1] = this->ComputePodY(vertex[1], normal, d, true);
    pod[2] = vertex[2];

    if (IsInsideImage(pod, 2)) {
      double intensity = InterpolateImageIntensity(pod);
      bgs.push_back(intensity);
    }
  }

  if (bgs.empty())
    return -1.0;
  else
    return Mean(bgs);
}



double Snake::ComputeForegroundMeanIntensity(unsigned index) const {
  const PointType &vertex = vertices_[index];
  const VectorType &normal = this->ComputeUnitTangentVector(index);

  assert(std::fabs(normal.GetNorm() - 1.0) < 1e-9);
  VectorType z_axis;
  z_axis[0] = 0.0;
  z_axis[1] = 0.0;
  z_axis[2] = 1.0;
  VectorType projection = z_axis - normal[2] * normal;
  VectorType long_axis;
  long_axis[0] = 1.0;
  long_axis[1] = 0.0;
  long_axis[2] = 0.0;
  VectorType short_axis;
  short_axis[0] = 0.0;
  short_axis[1] = 1.0;
  short_axis[2] = 0.0;
  double projection_length = projection.GetNorm();
  if (projection_length > kEpsilon) {
    long_axis = projection / projection_length;
    short_axis = itk::CrossProduct(long_axis, normal);
    short_axis.Normalize();
  }
  DataContainer bgs;
  const double angle_step = 2 * kPi / number_of_sectors_;
  for (int r = 0; r < radial_save_foreground_; r++) {
    for (int s = 0; s < number_of_sectors_; s++) {
      double angle = s * angle_step;
      VectorType v = static_cast<double>(r) * (std::cos(angle) * long_axis +
                                               std::sin(angle) * short_axis);
      v[2] *= z_spacing_;
      PointType p = vertex + v;
      if (IsInsideImage(p)) {
        double intensity = InterpolateImageIntensity(p);
        if (intensity > background_)
          bgs.push_back(intensity);
      }
    }
  }

  if (bgs.empty()) {
    return -1.0;  // return a negative value intentionally
  } else {
    return Mean(bgs);
  }
}

double Snake::ComputeForegroundMeanIntensity2d(unsigned index) const {
  const VectorType &normal = this->ComputeUnitTangentVector(index);
  PointType vertex = vertices_[index];
  DataContainer bgs;

  for (int d = 0; d < radial_save_foreground_; d++) {
    PointType pod;
    pod[0] = this->ComputePodX(vertex[0], normal, d, true);
    pod[1] = this->ComputePodY(vertex[1], normal, d, false);
    pod[2] = vertex[2];

    if (IsInsideImage(pod, 2)) {
      double intensity = InterpolateImageIntensity(pod);
      bgs.push_back(intensity);
    }

    pod[0] = this->ComputePodX(vertex[0], normal, d, false);
    pod[1] = this->ComputePodY(vertex[1], normal, d, true);
    pod[2] = vertex[2];

    if (IsInsideImage(pod, 2)) {
      double intensity = InterpolateImageIntensity(pod);
      bgs.push_back(intensity);
    }
  }

  if (bgs.empty())
    return -1.0;
  else
    return Mean(bgs);
}

InterpolatorOutputType Snake::InterpolateImageIntensity(PointType coords) const {
  const ImageType::SizeType &size =
      image_->GetLargestPossibleRegion().GetSize();

  // Make sure the coords is within 0.5 pixels of image edge in each direction
  if (coords[0] < -0.5)          coords[0] = -0.5;
  if (coords[0] > size[0] - 0.5) coords[0] = size[0] - 0.5;

  if (coords[1] < -0.5)          coords[1] = -0.5;
  if (coords[1] > size[1] - 0.5) coords[1] = size[1] - 0.5;

  if (coords[2] < - 0.5)         coords[2] = -0.5;
  if (coords[2] > size[2] - 0.5) coords[2] = size[2] - 0.5;

  return interpolator_->Evaluate(coords);
}

double Snake::ComputePodX(double x, const VectorType &tvec,
                          double dist, bool plus_root) const {
  if (std::abs(tvec[0]) < kEpsilon) {
    if (plus_root)
      return x + dist;
    else
      return x - dist;
  } else {
    double frac = tvec[1] / tvec[0];
    if (plus_root)
      return x + frac * dist / std::sqrt(1 + frac * frac);
    else
      return x - frac * dist / std::sqrt(1 + frac * frac);
  }
}


double Snake::ComputePodY(double y, const VectorType &tvec,
                          double dist, bool plus_root) const {
  if (std::abs(tvec[0]) < kEpsilon) {
    return y;
  } else {
    double frac = tvec[1] / tvec[0];
    if (plus_root)
      return y + dist / std::sqrt(1 + frac * frac);
    else
      return y - dist / std::sqrt(1 + frac * frac);
  }
}

bool Snake::CheckOrthogonality(const VectorType &vec1,
                               const VectorType &vec2) const {
  const double angle = 90;

  double dotprod = vec1 * vec2;
  double result = std::acos(dotprod) * 180.0 / kPi;

  if (std::abs(result-angle) < 1.0)
    return true;
  else
    return false;
}


void Snake::GetStartingRadialDirection(VectorType &direction,
                                       const VectorType &normal,
                                       const PointType &vertex) const {
  unsigned principal_index = this->GetPrincipalIndex(normal);
  // std::cout << "principal_index: " << principal_index << std::endl;
  PointType other;
  other[(principal_index + 1) % 3] = 0;
  other[(principal_index + 2) % 3] = 0;
  other[principal_index] = (normal[0]*vertex[0] + normal[1]*vertex[1] +
                            normal[2]*vertex[2]) / normal[principal_index];
  // std::cout << "other vertex: " << other << std::endl;
  // std::cout << "vertex: " << vertex << std::endl;
  direction = other - vertex;
  direction.Normalize();
}


unsigned Snake::GetPrincipalIndex(const VectorType &vec) const {
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
                               int d, int s) const {
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

bool Snake::IsInsideImage(const PointType &point, unsigned dim,
                          double padding) const {
  ImageType::SizeType size = image_->GetLargestPossibleRegion().GetSize();
  for (unsigned i = 0; i < dim; ++i) {
    if (point[i] < padding || point[i] > size[i] - padding)
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
void Snake::CheckBodyOverlap(const SnakeContainer &converged_snakes, const std::vector<std::vector<IndexPairContainer > >  &converged_snakes_grid) {
  if (!converged_) return;

  bool last_is_overlap = true;
  PointIterator overlap_start = vertices_.end();
  PointIterator overlap_end = vertices_.end();
  // PointType overlap_point;
  for (PointIterator it = vertices_.begin(); it != vertices_.end(); ++it) {
    bool overlap = Snake::VertexOverlap(*it, converged_snakes, converged_snakes_grid);
    if (overlap && !last_is_overlap) {
      overlap_start = it;
    } else if (!overlap && last_is_overlap && it > overlap_start) {
      overlap_end = it;
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
  std::cout << "snake id: " << this << std::endl;
  std::cout << "viable: " << viable_ << std::endl;
  std::cout << "open: " << std::boolalpha << open_
            << std::endl << std::noboolalpha;
  std::cout << "iteration: " << iterations_ << std::endl;
  std::cout << "length: " << length_ << std::endl;
  std::cout << "size: " << this->GetSize() << std::endl;
  std::cout << "spacing: " << spacing_ << std::endl;

  std::cout << "intensity: " << this->ComputeIntensity() << std::endl;
  std::cout << "local snr: " << this->ComputeSNR() << std::endl;

  std::cout << "fixed head: " << fixed_head_ << std::endl;
  std::cout << "fixed tail: " << fixed_tail_ << std::endl;

  const unsigned column_width = 15;
  std::cout << "#" << std::endl;
  for (unsigned j = 0; j != vertices_.size(); ++j) {
    std::cout << j << "\t";
    std::cout << std::setw(column_width) << this->GetX(j)
              << std::setw(column_width) << this->GetY(j)
              << std::setw(column_width) << this->GetZ(j)
              << std::setw(column_width)
              << InterpolateImageIntensity(this->GetPoint(j))
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
    else
      delete snake;

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
  else
    delete snake;
}

void Snake::EvolveWithTipFixed(SolverBank *solver, unsigned max_iter,
                               unsigned dim) {
  unsigned iter = 0;
  fixed_head_ = vertices_.front();
  fixed_tail_ = vertices_.back();

  while (viable_ && iter < max_iter) {
    if (!(iterations_ % check_period_)) {
      if (this->IsConverged())
        break;
    }

    this->IterateOnce(solver, dim);
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

bool Snake::ComputeLocalSNRAtIndex(unsigned index, int radial_near,
                                   int radial_far, double &local_snr) const {
  double foreground = InterpolateImageIntensity(this->GetPoint(index));
  double bg_mean(0.0), bg_std(0.0);
  bool local_bg_defined = this->ComputeLocalBackgroundMeanStd(
      index, radial_near, radial_far, bg_mean, bg_std);

  if (local_bg_defined) {
    if (bg_std < kEpsilon) {
      if (foreground > bg_mean)
        local_snr = kPlusInfinity;
      else
        local_snr = 0.0;
    } else {
      local_snr = (foreground - bg_mean) / bg_std;
    }
  }

  return local_bg_defined;
}

bool Snake::ComputeLocalBackgroundMeanStd(unsigned index, int radial_near,
                                          int radial_far, double &mean,
                                          double &std) const {
  DataContainer bgs;
  const VectorType normal = this->ComputeUnitTangentVector(index);
  PointType vertex = vertices_.at(index);
  VectorType z_axis;
  z_axis[0] = 0.0;
  z_axis[1] = 0.0;
  z_axis[2] = 1.0;
  VectorType projection = z_axis - normal[2] * normal;
  VectorType long_axis;
  long_axis[0] = 1.0;
  long_axis[1] = 0.0;
  long_axis[2] = 0.0;
  VectorType short_axis;
  short_axis[0] = 0.0;
  short_axis[1] = 1.0;
  short_axis[2] = 0.0;
  double projection_length = projection.GetNorm();
  if (projection_length > kEpsilon) {
    long_axis = projection / projection_length;
    short_axis = itk::CrossProduct(long_axis, normal);
    short_axis.Normalize();
  }

  const int number_of_sectors = 16;
  const double angle_step = 2 * kPi / number_of_sectors;
  for (int r = radial_near; r < radial_far; r++) {
    for (int s = 0; s < number_of_sectors; s++) {
      double angle = s * angle_step;
      VectorType v = static_cast<double>(r) * (std::cos(angle) * long_axis +
                                               std::sin(angle) * short_axis);
      v[2] *= z_spacing_;
      PointType p = vertex + v;
      if (IsInsideImage(p)) {
        bgs.push_back(InterpolateImageIntensity(p));
      }
    }
  }

  bool local_bg_defined = bgs.size() > number_of_sectors / 2;
  if (local_bg_defined) {
    mean = Mean(bgs);
    std = StandardDeviation(bgs, mean);
  }
  return local_bg_defined;
}

VectorType Snake::ComputeUnitTangentVector(unsigned index) const {
  VectorType tangent;

  if (index == 0) {
    if (head_tangent_[0]) {
      return head_tangent_;
    } else {
      tangent = vertices_.at(0) - vertices_.at(1);
    }
  } else if (index == vertices_.size() - 1) {
    if (tail_tangent_[0]) {
      return tail_tangent_;
    } else {
      tangent = vertices_.at(vertices_.size() - 2) -
                vertices_.at(vertices_.size() - 1);
    }
  } else {
    tangent = vertices_.at(index - 1) - vertices_.at(index + 1);
  }

  tangent.Normalize();
  return tangent;
}

void Snake::Trim(unsigned start, unsigned end) {
  vertices_.erase(vertices_.begin() + start,
                  vertices_.begin() + end);
}

void Snake::ExtendHead(const PointType &p) {
  vertices_.push_front(p);
}

void Snake::ExtendTail(const PointType &p) {
  vertices_.push_back(p);
}

void Snake::TrimAndInsert(unsigned start, unsigned end, const PointType &p) {
  if (start > end) {
    unsigned temp = start;
    start = end;
    end = temp;
  }
  vertices_.insert(vertices_.erase(vertices_.begin() + start,
                                   vertices_.begin() + end), p);
}

double Snake::ComputeIntensity() const {
  double intensity_sum = 0.0;
  for (unsigned i = 0; i < vertices_.size(); ++i) {
    intensity_sum += InterpolateImageIntensity(vertices_.at(i));
  }
  return intensity_sum / vertices_.size();
}

double Snake::ComputeSNR() const {
  double sum = 0.0;
  unsigned cnt = 0;
  for (unsigned i = 0; i < vertices_.size(); i++) {
    double snr = 0.0;
    bool bg_exist = this->ComputeLocalSNRAtIndex(
        i, radial_near_, radial_far_, snr);
    if (bg_exist) {
      sum += snr;
      cnt++;
    }
  }

  if (cnt)
    return sum / cnt;
  else
    return 0.0;
}

}  // namespace soax
