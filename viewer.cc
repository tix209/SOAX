/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file implements the visulization class for SOAX.
 */

#include "./viewer.h"
#include "QVTKWidget.h"
#include "itkImageToVTKImageFilter.h"
#include "vtkImageCast.h"
#include "vtkImagePlaneWidget.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"
#include "vtkVolume.h"
#include "vtkOrientationMarkerWidget.h"
#include "vtkCubeAxesActor.h"
#include "vtkCornerAnnotation.h"
#include "vtkVolumeProperty.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkVolumeRayCastMIPFunction.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkAxesActor.h"
#include "vtkCaptionActor2D.h"
#include "vtkTextProperty.h"
#include "vtkProperty.h"
#include "vtkPolyDataMapper.h"
#include "vtkOutlineSource.h"
#include "vtkDataSetMapper.h"
#include "vtkCellArray.h"
#include "vtkSphereSource.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"
#include "vtkTIFFWriter.h"
#include "vtkGL2PSExporter.h"
#include "vtkEventQtSlotConnect.h"
#include "vtkPointPicker.h"
#include "itkStatisticsImageFilter.h"
#include "./snake.h"
#include "./utility.h"

namespace soax {

double Viewer::kWhite[3] = {1.0, 1.0, 1.0};
double Viewer::kGray[3] = {0.8, 0.8, 0.8};
double Viewer::kRed[3] = {1.0, 0.0, 0.0};
double Viewer::kMagenta[3] = {1.0, 0.0, 1.0};
double Viewer::kYellow[3] = {1.0, 1.0, 0.0};
double Viewer::kGreen[3] = {0.0, 1.0, 0.0};
double Viewer::kCyan[3] = {0.0, 1.0, 1.0};
double Viewer::kBlue[3] = {0.0, 0.0, 1.0};

Viewer::Viewer():
    window_(0.0), level_(0.0), mip_min_intensity_(0.0),
    mip_max_intensity_(0.0), clip_span_(6.0), color_segment_step_(3),
    trimmed_actor_(NULL), selected_snake_(NULL), trimmed_snake_(NULL),
    trim_tip_index_(kBigNumber), trim_body_index1_(kBigNumber),
    trim_body_index2_(kBigNumber), is_trim_body_second_click_(false) {
  qvtk_ = new QVTKWidget;
  renderer_ = vtkSmartPointer<vtkRenderer>::New();
  qvtk_->GetRenderWindow()->AddRenderer(renderer_);
  camera_ = vtkSmartPointer<vtkCamera>::New();
  renderer_->SetActiveCamera(camera_);

  for (unsigned i = 0; i < kDimension; i++) {
    slice_planes_[i] = vtkImagePlaneWidget::New();
    slice_planes_[i]->SetResliceInterpolateToCubic();
    slice_planes_[i]->DisplayTextOn();
    slice_planes_[i]->SetInteractor(qvtk_->GetInteractor());
    slice_planes_[i]->PlaceWidget();
    slice_planes_[i]->SetSliceIndex(0);
    slice_planes_[i]->SetMarginSizeX(0);
    slice_planes_[i]->SetMarginSizeY(0);
    slice_planes_[i]->SetRightButtonAction(
        vtkImagePlaneWidget::VTK_SLICE_MOTION_ACTION);
    slice_planes_[i]->SetMiddleButtonAction(
        vtkImagePlaneWidget::VTK_WINDOW_LEVEL_ACTION);
  }

  volume_ = vtkSmartPointer<vtkVolume>::New();
  orientation_marker_ = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
  corner_text_ = vtkSmartPointer<vtkCornerAnnotation>::New();
  cube_axes_ = vtkSmartPointer<vtkCubeAxesActor>::New();
  bounding_box_ = vtkSmartPointer<vtkActor>::New();
  slot_connector_ = vtkSmartPointer<vtkEventQtSlotConnect>::New();
  clipped_actor_ = vtkSmartPointer<vtkActor>::New();
  picker_ = vtkSmartPointer<vtkPointPicker>::New();
  picker_->SetTolerance(0.01);
  qvtk_->GetInteractor()->SetPicker(picker_);
  snakes_actor_ = NULL;
  on_snake_sphere1_ = vtkActor::New();
  on_snake_sphere2_ = vtkActor::New();
  off_snake_sphere_ = vtkActor::New();
  inserted_point_.Fill(-1.0);

  snake_color_ = kMagenta;
  comparing_snakes1_color_ = kYellow;
  comparing_snakes2_color_ = kCyan;
  selected_snake_color_ = kCyan;
  snake_width_ = 3.0;
  comparing_snakes1_width_ = 8.0;
  comparing_snakes2_width_ = 14.0;
  snake_opacity_ = 1.0;
  comparing_snakes1_opacity_ = 0.5;
  comparing_snakes2_opacity_ = 0.3;
  junction_radius_ = 1.0;
  junction_color_ = kGreen;
  selected_junction_color_ = kBlue;
  sphere_color_ = kRed;
}

Viewer::~Viewer() {
  selected_junctions_.clear();
  selected_snakes_.clear();
  this->RemoveJunctions();
  this->RemoveSnakes();
  this->RemoveColorSegments();
  on_snake_sphere1_->Delete();
  on_snake_sphere2_->Delete();
  off_snake_sphere_->Delete();
  for (unsigned i = 0; i < kDimension; i++) {
    slice_planes_[i]->Delete();
  }
  delete qvtk_;
}

void Viewer::Reset() {
  window_ = level_ = mip_min_intensity_ = mip_max_intensity_ = 0.0;
  snake_filename_ = comparing_snake_filename1_ =
      comparing_snake_filename2_ = "";
  this->RemoveSnakes();
  this->RemoveJunctions();
  orientation_marker_->SetEnabled(false);
  corner_text_->ClearAllTexts();
  renderer_->RemoveActor(cube_axes_);
  renderer_->RemoveActor(bounding_box_);
  this->RemoveColorSegments();
  this->ResetTrimTip();
  this->ResetExtendTip();
  this->ResetTrimBody();
  renderer_->ResetCamera();
}

void Viewer::SetupImage(ImageType::Pointer image) {
  typedef itk::ImageToVTKImageFilter<ImageType>  ConnectorType;
  ConnectorType::Pointer connector = ConnectorType::New();
  connector->SetInput(image);
  connector->Update();

  typedef itk::StatisticsImageFilter<ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->Update();
  double min_intensity = filter->GetMinimum();
  double max_intensity = filter->GetMaximum();
  window_ = max_intensity - min_intensity;
  level_ = min_intensity + window_ / 2;
  mip_min_intensity_ = min_intensity +
      0.05 * (max_intensity - min_intensity);
  mip_max_intensity_ = max_intensity;

  this->SetupSlicePlanes(connector->GetOutput());
  this->SetupMIPRendering(connector->GetOutput());
  this->SetupOrientationMarker();
  this->SetupUpperLeftCornerText(min_intensity, max_intensity);
  this->SetupBoundingBox(volume_);
  this->SetupCubeAxes(connector->GetOutput());

  this->UpdateJunctionRadius(image);
}

void Viewer::UpdateJunctionRadius(ImageType::Pointer image) {
  ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
  double diag_length = 0.0;
  for (unsigned i = 0; i < kDimension; ++i) {
    diag_length += size[i] * size[i];
  }
  diag_length = std::sqrt(diag_length);
  if (diag_length < 100.0)
    junction_radius_ = 1.5;
  else if (diag_length < 200.0)
    junction_radius_ = 2;
  else if (diag_length < 500.0)
    junction_radius_ = 2.5;
  else
    junction_radius_ = 3.5;
}

void Viewer::SetupSlicePlanes(vtkImageData *data) {
  for (unsigned i = 0; i < kDimension; ++i) {
    slice_planes_[i]->SetInputData(data);
    slice_planes_[i]->SetPlaneOrientation(i);
    slice_planes_[i]->UpdatePlacement();
    slice_planes_[i]->SetWindowLevel(window_, level_);
  }
}

void Viewer::ToggleSlicePlanes(bool state) {
  if (state) {
    for (unsigned i = 0; i < kDimension; ++i) {
      slice_planes_[i]->On();
    }
    renderer_->ResetCamera();
  } else {
    for (unsigned i = 0; i < kDimension; ++i) {
      slice_planes_[i]->Off();
    }
  }
  this->Render();
}

void Viewer::SetupMIPRendering(vtkImageData *data) {
  vtkSmartPointer<vtkPiecewiseFunction> opacity_function =
      vtkSmartPointer<vtkPiecewiseFunction>::New();

  opacity_function->AddPoint(mip_min_intensity_, 0.0);
  opacity_function->AddPoint(mip_max_intensity_, 1.0);

  vtkSmartPointer<vtkColorTransferFunction> color_function =
      vtkSmartPointer<vtkColorTransferFunction>::New();
  color_function->SetColorSpaceToRGB();
  // color is all white
  color_function->AddRGBPoint(0, 1, 1, 1);
  color_function->AddRGBPoint(255, 1, 1, 1);

  vtkSmartPointer<vtkVolumeProperty> mip_volume_property =
      vtkSmartPointer<vtkVolumeProperty>::New();
  mip_volume_property->SetScalarOpacity(opacity_function);
  mip_volume_property->SetColor(color_function);
  mip_volume_property->SetInterpolationTypeToLinear();
  volume_->SetProperty(mip_volume_property);

  vtkSmartPointer<vtkVolumeRayCastMapper> mapper =
      vtkSmartPointer<vtkVolumeRayCastMapper>::New();
  mapper->SetInputData(data);
  vtkSmartPointer<vtkVolumeRayCastMIPFunction> mip_function =
      vtkSmartPointer<vtkVolumeRayCastMIPFunction>::New();
  mapper->SetVolumeRayCastFunction(mip_function);
  volume_->SetMapper(mapper);
}

void Viewer::ToggleMIPRendering(bool state) {
  if (state) {
    renderer_->AddViewProp(volume_);
    renderer_->ResetCamera();
  } else {
    renderer_->RemoveViewProp(volume_);
  }
  this->Render();
}

void Viewer::UpdateWindowLevel(double window, double level) {
  window_ = window;
  level_ = level;
  for (unsigned i = 0; i < kDimension; ++i) {
    slice_planes_[i]->SetWindowLevel(window, level);
  }
}

void Viewer::UpdateMIPIntensityRange(double min, double max) {
  mip_min_intensity_ = min;
  mip_max_intensity_ = max;
  volume_->GetProperty()->GetScalarOpacity()->RemoveAllPoints();
  volume_->GetProperty()->GetScalarOpacity()->AddPoint(min, 0.0);
  volume_->GetProperty()->GetScalarOpacity()->AddPoint(max, 1.0);
}


void Viewer::SetupOrientationMarker() {
  vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
  const int font_size = 14;
  axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->
      SetFontSize(font_size);
  axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->
      SetFontSize(font_size);
  axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->
      SetFontSize(font_size);
  orientation_marker_->SetOrientationMarker(axes);
  orientation_marker_->SetInteractor(qvtk_->GetInteractor());
  orientation_marker_->SetViewport(0, 0, 0.3, 0.3);
}

void Viewer::ToggleOrientationMarker(bool state) {
  orientation_marker_->SetEnabled(state);
  if (state)
    orientation_marker_->InteractiveOff();
  this->Render();
}

void Viewer::SetupUpperLeftCornerText(unsigned min_intensity,
                                      unsigned max_intensity) {
  std::ostringstream buffer;
  buffer << "Intensity: [" << min_intensity << ", "
         << max_intensity << "]";
  corner_text_->SetText(2, buffer.str().c_str());
}

void Viewer::SetupUpperRightCornerText() {
  std::string info;
  if (!snake_filename_.empty())
    info += snake_filename_ + " (Magenta)\n";
  if (!comparing_snake_filename1_.empty())
    info += comparing_snake_filename1_ + " (Yellow)\n";
  if (!comparing_snake_filename2_.empty())
    info += comparing_snake_filename2_ + " (Cyan)";
  corner_text_->SetText(3, info.c_str());
}

void Viewer::ToggleCornerText(bool state) {
  if (state)
    renderer_->AddViewProp(corner_text_);
  else
    renderer_->RemoveViewProp(corner_text_);
  this->Render();
}

void Viewer::SetupBoundingBox(vtkSmartPointer<vtkVolume> volume) {
  vtkSmartPointer<vtkOutlineSource> outline =
      vtkSmartPointer<vtkOutlineSource>::New();
  vtkSmartPointer<vtkPolyDataMapper> outline_mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  outline_mapper->SetInputConnection(outline->GetOutputPort());
  outline->SetBounds(volume->GetBounds());
  bounding_box_->PickableOff();
  bounding_box_->DragableOff();
  bounding_box_->SetMapper(outline_mapper);
  bounding_box_->GetProperty()->SetLineWidth(1.0);
  bounding_box_->GetProperty()->SetColor(kRed);
  bounding_box_->GetProperty()->SetAmbient(1.0);
  bounding_box_->GetProperty()->SetDiffuse(0.0);
}

void Viewer::ToggleBoundingBox(bool state) {
  if (state)
    renderer_->AddActor(bounding_box_);
  else
    renderer_->RemoveActor(bounding_box_);
  this->Render();
}

void Viewer::SetupCubeAxes(vtkSmartPointer<vtkImageData> image) {
  cube_axes_->SetBounds(image->GetBounds());
  cube_axes_->SetCamera(camera_);
  cube_axes_->GetTitleTextProperty(0)->SetColor(kRed);
  cube_axes_->GetLabelTextProperty(0)->SetColor(kRed);
  cube_axes_->GetTitleTextProperty(1)->SetColor(kGreen);
  cube_axes_->GetLabelTextProperty(1)->SetColor(kGreen);
  cube_axes_->GetTitleTextProperty(2)->SetColor(kBlue);
  cube_axes_->GetLabelTextProperty(2)->SetColor(kBlue);
  vtkSmartPointer<vtkProperty> axes_property =
      vtkSmartPointer<vtkProperty>::New();
  axes_property->SetLineWidth(2.0);
  axes_property->SetColor(kRed);
  cube_axes_->SetXAxesLinesProperty(axes_property);
  axes_property->SetColor(kGreen);
  cube_axes_->SetYAxesLinesProperty(axes_property);
  axes_property->SetColor(kBlue);
  cube_axes_->SetZAxesLinesProperty(axes_property);
}

void Viewer::ToggleCubeAxes(bool state) {
  if (state)
    renderer_->AddActor(cube_axes_);
  else
    renderer_->RemoveActor(cube_axes_);
  this->Render();
}

void Viewer::SetupSnakesAsOneActor(const SnakeContainer &snakes) {
  if (snakes.empty()) return;
  vtkPolyData *curve = this->MakePolyDataForMultipleSnakes(snakes);
  vtkDataSetMapper *mapper = vtkDataSetMapper::New();
  mapper->SetInputData(curve);
  snakes_actor_ = vtkActor::New();
  snakes_actor_->SetMapper(mapper);
  mapper->Update();
  mapper->Delete();
  curve->Delete();
  this->SetupEvolvingActorProperty(snakes_actor_);
  renderer_->AddActor(snakes_actor_);
}

void Viewer::SetupSnakes(const SnakeContainer &snakes, unsigned category) {
  if (snakes.empty()) return;

  for (SnakeConstIterator it = snakes.begin(); it != snakes.end(); ++it) {
    this->SetupSnake(*it, category);
  }
}

void Viewer::SetupSnake(Snake *snake, unsigned category) {
  // Remove the old snake from the scene and actor-snake map
  SnakeActorMap::iterator it = snake_actor_map_.find(snake);
  if (it != snake_actor_map_.end()) {
    renderer_->RemoveActor(it->second);
    actor_snake_map_.erase(it->second);
    it->second->Delete();
  }

  vtkActor *actor = this->ActSnake(snake);
  actor_snake_map_[actor] = snake;
  snake_actor_map_[snake] = actor;

  switch (category) {
    case 0:
      this->SetupEvolvingActorProperty(actor);
      break;
    case 1:
      this->SetupComparingActorProperty(actor);
      break;
    case 2:
      this->SetupAnotherComparingActorProperty(actor);
      break;
    default:
      std::cerr << "SetupSnake: unknown snake category!" << std::endl;
  }
  renderer_->AddActor(actor);
}

vtkActor * Viewer::ActSnake(Snake *snake) {
  return this->ActSnakeSegments(snake, 0, snake->GetSize());
}

vtkActor * Viewer::ActSnakeSegments(Snake *snake,
                                    unsigned start, unsigned end) {
  vtkDataSetMapper *mapper = vtkDataSetMapper::New();
  vtkPolyData *curve = this->MakePolyData(snake, start, end);
  mapper->SetInputData(curve);
  vtkActor *actor = vtkActor::New();
  actor->SetMapper(mapper);
  mapper->Update();
  mapper->Delete();
  curve->Delete();
  return actor;
}


vtkPolyData * Viewer::MakePolyData(Snake *snake,
                                   unsigned start, unsigned end) {
  vtkPolyData *curve = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *cells = vtkCellArray::New();

  vtkIdType cell_index[2];
  vtkFloatingPointType coordinates[3];

  for (unsigned i = start, j = 0; i < end; ++i, ++j) {
    coordinates[0] = snake->GetX(i);
    coordinates[1] = snake->GetY(i);
    coordinates[2] = snake->GetZ(i);

    points->InsertPoint(j, coordinates);

    // last point
    if (i == end - 1) {
      if (!snake->open()) {
        cell_index[0] = end-1-start;
        cell_index[1] = 0;
        cells->InsertNextCell(2, cell_index);
      }
    } else {
      cell_index[0] = j;
      cell_index[1] = j+1;
      cells->InsertNextCell(2, cell_index);
    }
  }

  curve->SetPoints(points);
  curve->SetLines(cells);
  points->Delete();
  cells->Delete();
  return curve;
}

vtkPolyData * Viewer::MakePolyDataForMultipleSnakes(
    const SnakeContainer &snakes) {
  vtkPolyData *curve = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *cells = vtkCellArray::New();

  unsigned int index = 0;
  vtkIdType cell_index[2];
  vtkFloatingPointType coordinates[3];

  for (SnakeContainer::const_iterator s_iter = snakes.begin();
       s_iter != snakes.end(); ++s_iter) {
    for (unsigned i = 0; i < (*s_iter)->GetSize(); ++i) {
      coordinates[0] = (*s_iter)->GetX(i);
      coordinates[1] = (*s_iter)->GetY(i);
      coordinates[2] = (*s_iter)->GetZ(i);
      points->InsertPoint(index, coordinates);
      if (i != (*s_iter)->GetSize() - 1) {
        cell_index[0] = index;
        cell_index[1] = index + 1;
        cells->InsertNextCell(2, cell_index);
      }
      index++;
    }
  }
  curve->SetPoints(points);
  curve->SetLines(cells);
  points->Delete();
  cells->Delete();

  return curve;
}

void Viewer::SetupEvolvingActorProperty(vtkActor *actor) {
  actor->GetProperty()->SetInterpolationToPhong();
  actor->GetProperty()->SetOpacity(snake_opacity_);
  actor->GetProperty()->SetAmbient(0.2);
  actor->GetProperty()->SetDiffuse(0.7);
  actor->GetProperty()->SetSpecular(0.6);
  actor->GetProperty()->SetSpecularPower(50);
  actor->GetProperty()->SetColor(snake_color_);
  actor->GetProperty()->SetLineWidth(snake_width_);
}

void Viewer::SetupComparingActorProperty(vtkActor *actor) {
  actor->GetProperty()->SetInterpolationToPhong();
  actor->GetProperty()->SetOpacity(comparing_snakes1_opacity_);
  actor->GetProperty()->SetAmbient(0.2);
  actor->GetProperty()->SetDiffuse(0.7);
  actor->GetProperty()->SetSpecular(0.6);
  actor->GetProperty()->SetSpecularPower(50);
  actor->GetProperty()->SetColor(comparing_snakes1_color_);
  actor->GetProperty()->SetLineWidth(comparing_snakes1_width_);
}

void Viewer::SetupAnotherComparingActorProperty(vtkActor *actor) {
  actor->GetProperty()->SetInterpolationToPhong();
  actor->GetProperty()->SetOpacity(comparing_snakes2_opacity_);
  actor->GetProperty()->SetAmbient(0.2);
  actor->GetProperty()->SetDiffuse(0.7);
  actor->GetProperty()->SetSpecular(0.6);
  actor->GetProperty()->SetSpecularPower(50);
  actor->GetProperty()->SetColor(comparing_snakes2_color_);
  actor->GetProperty()->SetLineWidth(comparing_snakes2_width_);
}

void Viewer::ChangeSnakeColor(Snake *s, double *color) {
  snake_actor_map_[s]->GetProperty()->SetColor(color);
}

void Viewer::ToggleSnakes(bool state) {
  if (state) {
    if (snakes_actor_)
      renderer_->AddActor(snakes_actor_);
    for (ActorSnakeMap::iterator it = actor_snake_map_.begin();
         it != actor_snake_map_.end(); ++it) {
      renderer_->AddActor(it->first);
    }
  } else {
    renderer_->RemoveActor(snakes_actor_);
    for (ActorSnakeMap::iterator it = actor_snake_map_.begin();
         it != actor_snake_map_.end(); ++it) {
      renderer_->RemoveActor(it->first);
    }
  }
  this->Render();
}

void Viewer::RemoveSnakes() {
  for (ActorSnakeMap::iterator it = actor_snake_map_.begin();
       it != actor_snake_map_.end(); ++it) {
    renderer_->RemoveActor(it->first);
    it->first->Delete();
  }
  snake_actor_map_.clear();
  actor_snake_map_.clear();
  selected_snakes_.clear();
  if (snakes_actor_) {
    renderer_->RemoveActor(snakes_actor_);
    snakes_actor_->Delete();
    snakes_actor_ = NULL;
  }
}

void Viewer::RemoveSnake(Snake *snake) {
  SnakeActorMap::iterator it = snake_actor_map_.find(snake);
  if (it == snake_actor_map_.end()) {
    return;
  }

  vtkActor *actor = it->second;
  renderer_->RemoveActor(actor);
  snake_actor_map_.erase(snake);
  actor_snake_map_.erase(actor);
  actor->Delete();
}

void Viewer::SetupJunctions(const PointContainer &points) {
  if (points.empty()) return;
  for (PointConstIterator it = points.begin(); it != points.end(); ++it) {
    vtkActor *s = vtkActor::New();
    this->SetupSphere(*it, s, junction_color_);
    actor_junctions_[s] = *it;
  }
}

void Viewer::SetupSphere(const PointType &point, vtkActor *sphere,
                         double *color) {
  vtkSphereSource *source = vtkSphereSource::New();
  source->SetCenter(point[0], point[1], point[2]);
  source->SetRadius(junction_radius_);
  source->Update();
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
  mapper->SetInputConnection(source->GetOutputPort());
  sphere->SetMapper(mapper);
  sphere->GetProperty()->SetColor(color);
  source->Delete();
  mapper->Delete();
  renderer_->AddActor(sphere);
}

void Viewer::ToggleJunctions(bool state) {
  if (state) {
    for (ActorPointMap::const_iterator it = actor_junctions_.begin();
         it != actor_junctions_.end(); ++it) {
      renderer_->AddActor(it->first);
    }
  } else {
    for (ActorPointMap::const_iterator it = actor_junctions_.begin();
         it != actor_junctions_.end(); ++it) {
      renderer_->RemoveActor(it->first);
    }
  }
  this->Render();
}

void Viewer::RemoveJunctions() {
  for (ActorPointMap::iterator it = actor_junctions_.begin();
       it != actor_junctions_.end(); ++it) {
    renderer_->RemoveActor(it->first);
    it->first->Delete();
  }
  actor_junctions_.clear();
  selected_junctions_.clear();
}

void Viewer::ToggleClipSnakes(bool state) {
  if (state) {
    for (unsigned i = 0; i < kDimension; ++i) {
      slot_connector_->Connect(
          slice_planes_[i], vtkCommand::EndInteractionEvent,
          this, SLOT(SetupClippedSnakes(vtkObject *)));
    }
  } else {
    for (unsigned i = 0; i < kDimension; ++i) {
      slot_connector_->Disconnect(
          slice_planes_[i], vtkCommand::EndInteractionEvent,
          this, SLOT(SetupClippedSnakes(vtkObject *)));
    }
    renderer_->RemoveActor(clipped_actor_);
  }
  this->Render();
}

void Viewer::SetupClippedSnakes(vtkObject *obj) {
  vtkImagePlaneWidget *plane = vtkImagePlaneWidget::SafeDownCast(obj);
  unsigned axis = plane->GetPlaneOrientation();
  double position = slice_planes_[axis]->GetSlicePosition();
  vtkPolyData *data = this->MakeClippedPolyData(axis, position);
  vtkDataSetMapper *mapper = vtkDataSetMapper::New();
  mapper->SetInputData(data);
  clipped_actor_->SetMapper(mapper);
  mapper->Update();
  data->Delete();
  this->SetupEvolvingActorProperty(clipped_actor_);
  renderer_->AddActor(clipped_actor_);
}

vtkPolyData * Viewer::MakeClippedPolyData(unsigned axis, double position) {
  vtkPolyData *curve = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *cells = vtkCellArray::New();

  unsigned int index = 0;
  vtkIdType cell_index[2];
  vtkFloatingPointType coordinates[3];

  for (SnakeActorMap::const_iterator it = snake_actor_map_.begin();
       it != snake_actor_map_.end(); ++it) {
    for (unsigned i = 0; i < it->first->GetSize(); ++i) {
      coordinates[0] = it->first->GetX(i);
      coordinates[1] = it->first->GetY(i);
      coordinates[2] = it->first->GetZ(i);

      if (coordinates[axis] < position - clip_span_ ||
          coordinates[axis] > position + clip_span_)
        continue;

      points->InsertPoint(index, coordinates);

      if (i != it->first->GetSize() - 1) {
        double next_coordinate = it->first->GetPoint(i+1)[axis];
        if (next_coordinate >= position - clip_span_ &&
            next_coordinate <= position + clip_span_) {
          cell_index[0] = index;
          cell_index[1] = index + 1;
          cells->InsertNextCell(2, cell_index);
        }
      }
      index++;
    }
  }
  curve->SetPoints(points);
  curve->SetLines(cells);
  points->Delete();
  cells->Delete();

  return curve;
}

void Viewer::ColorByAzimuthalAngle(bool state) {
  this->ColorSnakes(state, true);
}

void Viewer::ColorByPolarAngle(bool state) {
  this->ColorSnakes(state, false);
}

void Viewer::ColorSnakes(bool state, bool azimuthal) {
  if (state) {
    this->SetupColorSegments(azimuthal);
  } else {
    for (unsigned i = 0; i < color_segments_.size(); ++i) {
      renderer_->RemoveActor(color_segments_[i]);
    }
  }
  this->Render();
}

void Viewer::SetupColorSegments(bool azimuthal) {
  if (actor_snake_map_.empty()) return;
  for (ActorSnakeMap::const_iterator it = actor_snake_map_.begin();
       it != actor_snake_map_.end(); ++it) {
    unsigned step = color_segment_step_ > it->second->GetSize() - 1 ?
        it->second->GetSize() - 1 : color_segment_step_;
    unsigned i = 0;
    while (i < it->second->GetSize() - step) {
      vtkActor *actor = this->ActSnakeSegments(it->second, i,
                                               i + step + 1);
      this->SetupEvolvingActorProperty(actor);
      double color[3];
      VectorType vector = it->second->GetPoint(i) -
          it->second->GetPoint(i + step);
      this->ComputeColor(vector, azimuthal, color);
      actor->GetProperty()->SetColor(color);
      color_segments_.push_back(actor);
      renderer_->AddActor(actor);
      i += step;
    }
    if (i < it->second->GetSize() - 1) {
      vtkActor *actor = this->ActSnakeSegments(it->second, i,
                                               it->second->GetSize());
      this->SetupEvolvingActorProperty(actor);
      double color[3];
      VectorType vector = it->second->GetPoint(i) -
          it->second->GetTail();
      this->ComputeColor(vector, azimuthal, color);
      actor->GetProperty()->SetColor(color);
      color_segments_.push_back(actor);
      renderer_->AddActor(actor);
    }
  }
}

void Viewer::ComputeColor(const VectorType &vector,
                          bool azimuthal, double *color) {
  double theta(0.0), phi(0.0);
  this->ComputeThetaPhi(vector, &theta, &phi);
  double hue = 0.0;
  if (azimuthal) {
    hue = phi * 2 + 180;
  } else {
    hue = theta * 2;
  }
  this->ComputeRGBFromHue(hue, color);
}

void Viewer::ComputeThetaPhi(const VectorType &vector,
                             double *theta,
                             double *phi) {
  // phi is (-pi/2, +pi/2]
  // theta is [0, pi)
  VectorType v = vector;
  const double r = v.GetNorm();
  if (std::abs(v[0]) < kEpsilon && std::abs(v[1]) < kEpsilon) {
    // x = y = 0
    *phi = 0;
    *theta = 0;
  } else if (std::abs(v[0]) < kEpsilon) {
    // x = 0, y != 0
    if (v[1] < -kEpsilon)
      v = -v;

    *phi = 90;
    *theta = std::acos(v[2]/r) * 180 / kPi;
  } else {
    // x != 0
    if (v[0] < -kEpsilon)
      v = -v;

    *theta = std::acos(v[2]/r) * 180 / kPi;
    *phi = std::atan(v[1] / v[0]) * 180 / kPi;
  }
}

void Viewer::ComputeRGBFromHue(double hue, double *color) {
  double intensity = 0.5;
  if (hue < 120) {
    color[0] = intensity * (1 + std::cos(hue*kPi/180.0) /
                       std::cos((60-hue)*kPi/180));
    color[1] = 3*intensity - color[0];
    color[2] = 0;
  } else if (hue < 240) {
    double hue1 = hue - 120;
    color[0] = 0;
    color[1] = intensity * (1 + std::cos(hue1*kPi/180.0) /
                         std::cos((60-hue1)*kPi/180));
    color[2] = 3*intensity - color[1];
  } else {
    double hue1 = hue - 240;
    color[1] = 0;
    color[2] = intensity * (1 + std::cos(hue1*kPi/180.0) /
                        std::cos((60-hue1)*kPi/180));
    color[0] = 3*intensity - color[2];
  }
}

void Viewer::RemoveColorSegments() {
  for (ActorContainer::iterator it = color_segments_.begin();
       it != color_segments_.end(); ++it) {
    renderer_->RemoveActor(*it);
    (*it)->Delete();
  }
  color_segments_.clear();
}

void Viewer::ToggleNone(bool state) {
  if (state) {
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::LeftButtonPressEvent,
                             this, SLOT(SelectSnakeForView()));
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::RightButtonPressEvent,
                             this, SLOT(DeselectSnakeForView()));
  } else {
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::LeftButtonPressEvent,
                                this, SLOT(SelectSnakeForView()));
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::RightButtonPressEvent,
                                this, SLOT(DeselectSnakeForView()));
  }
}

void Viewer::SelectSnakeForView() {
  picker_->Pick(qvtk_->GetInteractor()->GetEventPosition()[0],
                qvtk_->GetInteractor()->GetEventPosition()[1],
                0, renderer_);
  vtkActor *actor = picker_->GetActor();
  if (actor) {
    ActorSnakeMap::const_iterator it = actor_snake_map_.find(actor);
    if (it != actor_snake_map_.end()) {
      actor->GetProperty()->SetColor(selected_snake_color_);
      selected_snake_ = it->second;
      it->second->PrintSelf();
    }
  }
}

void Viewer::DeselectSnakeForView() {
  picker_->Pick(qvtk_->GetInteractor()->GetEventPosition()[0],
                qvtk_->GetInteractor()->GetEventPosition()[1],
                0, renderer_);
  vtkActor *actor = picker_->GetActor();
  if (actor) {
    ActorSnakeMap::const_iterator it = actor_snake_map_.find(actor);
    if (it != actor_snake_map_.end()) {
      actor->GetProperty()->SetColor(snake_color_);
      selected_snake_ = NULL;
    }
  }
}

void Viewer::ToggleDeleteSnake(bool state) {
  if (state) {
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::LeftButtonPressEvent,
                             this, SLOT(SelectSnakeForDeletion()));
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::RightButtonPressEvent,
                             this, SLOT(DeselectSnakeForDeletion()));
  } else {
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::LeftButtonPressEvent,
                                this, SLOT(SelectSnakeForDeletion()));
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::RightButtonPressEvent,
                                this, SLOT(DeselectSnakeForDeletion()));
  }
}

void Viewer::SelectSnakeForDeletion() {
  picker_->Pick(qvtk_->GetInteractor()->GetEventPosition()[0],
                qvtk_->GetInteractor()->GetEventPosition()[1],
                0, renderer_);
  vtkActor *actor = picker_->GetActor();
  if (actor) {
    ActorSnakeMap::const_iterator it = actor_snake_map_.find(actor);
    if (it != actor_snake_map_.end()) {
      actor->GetProperty()->SetColor(selected_snake_color_);
      selected_snakes_.insert(it->second);
      it->second->PrintSelf();
    }
  }
}

void Viewer::DeselectSnakeForDeletion() {
  picker_->Pick(qvtk_->GetInteractor()->GetEventPosition()[0],
                qvtk_->GetInteractor()->GetEventPosition()[1],
                0, renderer_);
  vtkActor *actor = picker_->GetActor();
  if (actor) {
    ActorSnakeMap::const_iterator it = actor_snake_map_.find(actor);
    if (it != actor_snake_map_.end()) {
      actor->GetProperty()->SetColor(snake_color_);
      selected_snakes_.erase(it->second);
    }
  }
}

void Viewer::RemoveSelectedSnakes() {
  if (selected_snakes_.empty()) return;
  for (SnakeSet::const_iterator it = selected_snakes_.begin();
       it != selected_snakes_.end(); ++it) {
    SnakeActorMap::const_iterator snake_actor_it =
        snake_actor_map_.find(*it);
    if (snake_actor_it != snake_actor_map_.end()) {
      renderer_->RemoveActor(snake_actor_it->second);
      actor_snake_map_.erase(snake_actor_it->second);
      snake_actor_map_.erase(*it);
      snake_actor_it->second->Delete();
    }
  }
  selected_snakes_.clear();
}

void Viewer::ToggleTrimTip(bool state) {
  if (state) {
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::LeftButtonPressEvent,
                             this, SLOT(SelectVertex()));
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::RightButtonPressEvent,
                             this, SLOT(DeselectVertex()));
  } else {
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::LeftButtonPressEvent,
                                this, SLOT(SelectVertex()));
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::RightButtonPressEvent,
                                this, SLOT(DeselectVertex()));
    this->ResetTrimTip();
  }
  this->Render();
}

void Viewer::ResetTrimTip() {
  if (trimmed_actor_)
    trimmed_actor_->GetProperty()->SetColor(snake_color_);
  trimmed_actor_ = NULL;
  trimmed_snake_ = NULL;
  trim_tip_index_ = kBigNumber;
  renderer_->RemoveActor(on_snake_sphere1_);
}

void Viewer::SelectVertex() {
  picker_->Pick(qvtk_->GetInteractor()->GetEventPosition()[0],
                qvtk_->GetInteractor()->GetEventPosition()[1],
                0, renderer_);
  vtkActor *actor = picker_->GetActor();
  if (actor) {
    ActorSnakeMap::const_iterator it = actor_snake_map_.find(actor);
    if (it != actor_snake_map_.end()) {
      unsigned index = static_cast<unsigned>(picker_->GetPointId());
      this->SetupSphere(it->second->GetPoint(index), on_snake_sphere1_,
                        sphere_color_);
      it->first->GetProperty()->SetColor(selected_snake_color_);
      if (trimmed_actor_ && trimmed_actor_ != it->first) {
        // reset previously selected snake
        trimmed_actor_->GetProperty()->SetColor(snake_color_);
      }
      trimmed_actor_ = it->first;
      trimmed_snake_ = it->second;
      trim_tip_index_ = index;
      // trimmed_snake_->PrintSelf();
    }
  }
}

void Viewer::DeselectVertex() {
  picker_->Pick(qvtk_->GetInteractor()->GetEventPosition()[0],
                qvtk_->GetInteractor()->GetEventPosition()[1],
                0, renderer_);
  vtkActor *actor = picker_->GetActor();
  if (actor == on_snake_sphere1_) {
    renderer_->RemoveActor(on_snake_sphere1_);
    trimmed_actor_->GetProperty()->SetColor(snake_color_);
    trimmed_actor_ = NULL;
    trimmed_snake_ = NULL;
    trim_tip_index_ = kBigNumber;
  }
}

void Viewer::TrimTip() {
  if (!trimmed_snake_) return;
  unsigned mid = trimmed_snake_->GetSize()/2;
  if (trim_tip_index_ < mid) {
    trimmed_snake_->Trim(0, trim_tip_index_);
  } else {
    trimmed_snake_->Trim(trim_tip_index_, trimmed_snake_->GetSize());
  }
  this->SetupSnake(trimmed_snake_, 0);
  trim_tip_index_ = kBigNumber;
  trimmed_actor_ = NULL;
  renderer_->RemoveActor(on_snake_sphere1_);
}

void Viewer::ToggleExtendTip(bool state) {
  if (state) {
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::LeftButtonPressEvent,
                             this, SLOT(SelectVertex()));
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::RightButtonPressEvent,
                             this, SLOT(DeselectVertex()));

    for (unsigned i = 0; i < kDimension; ++i) {
      slot_connector_->Connect(slice_planes_[i],
                               vtkCommand::EndInteractionEvent,
                               this,
                               SLOT(SelectExtendVertex(vtkObject *)));
    }
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::RightButtonPressEvent,
                             this,
                             SLOT(DeselectExtendVertex(vtkObject *)));
  } else {
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::LeftButtonPressEvent,
                                this, SLOT(SelectVertex()));
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::RightButtonPressEvent,
                                this, SLOT(DeselectVertex()));
    for (unsigned i = 0; i < kDimension; ++i) {
      slot_connector_->Disconnect(slice_planes_[i],
                                  vtkCommand::EndInteractionEvent,
                                  this,
                                  SLOT(SelectExtendVertex(vtkObject *)));
    }
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::RightButtonPressEvent,
                                this, SLOT(DeselectVertex(vtkObject *)));
    this->ResetExtendTip();
  }
  this->Render();
}

void Viewer::ResetExtendTip() {
  if (trimmed_actor_)
    trimmed_actor_->GetProperty()->SetColor(snake_color_);
  trimmed_actor_ = NULL;
  trimmed_snake_ = NULL;
  renderer_->RemoveActor(on_snake_sphere1_);
  renderer_->RemoveActor(off_snake_sphere_);
  inserted_point_.Fill(-1.0);
}

void Viewer::SelectExtendVertex(vtkObject *obj) {
  vtkImagePlaneWidget *p = vtkImagePlaneWidget::SafeDownCast(obj);
  double point[3];
  p->GetCurrentCursorPosition(point);

  if (point[0] >= 0 && trim_tip_index_ != kBigNumber) {
    inserted_point_[0] = point[0];
    inserted_point_[1] = point[1];
    inserted_point_[2] = point[2];

    this->SetupSphere(inserted_point_, off_snake_sphere_, sphere_color_);
  }
}

void Viewer::DeselectExtendVertex(vtkObject *obj) {
  this->DeselectInsertedVertex(obj);
}

void Viewer::ExtendTip() {
  if (!trimmed_snake_ || inserted_point_[0] < 0) return;
  unsigned mid = trimmed_snake_->GetSize()/2;
  if (trim_tip_index_ < mid) {
    trimmed_snake_->ExtendHead(inserted_point_);
  } else {
    trimmed_snake_->ExtendTail(inserted_point_);
  }
  trimmed_snake_->Resample();
  this->SetupSnake(trimmed_snake_, 0);
  inserted_point_.Fill(-1.0);
  trim_tip_index_ = kBigNumber;
  trimmed_actor_ = NULL;
  renderer_->RemoveActor(on_snake_sphere1_);
  renderer_->RemoveActor(off_snake_sphere_);
}

void Viewer::ToggleTrimBody(bool state) {
  if (state) {
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::LeftButtonPressEvent,
                             this, SLOT(SelectBodyVertex()));
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::RightButtonPressEvent,
                             this, SLOT(DeselectBodyVertex()));

    for (unsigned i = 0; i < kDimension; ++i) {
      slot_connector_->Connect(slice_planes_[i],
                               vtkCommand::EndInteractionEvent,
                               this,
                               SLOT(SelectInsertedVertex(vtkObject *)));
    }
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::RightButtonPressEvent,
                             this,
                             SLOT(DeselectInsertedVertex(vtkObject *)));
  } else {
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::LeftButtonPressEvent,
                                this, SLOT(SelectBodyVertex()));
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::RightButtonPressEvent,
                                this, SLOT(DeselectBodyVertex()));
    for (unsigned i = 0; i < kDimension; ++i) {
      slot_connector_->Disconnect(slice_planes_[i],
                                  vtkCommand::EndInteractionEvent,
                                  this,
                                  SLOT(SelectInsertedVertex(vtkObject *)));
    }
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::RightButtonPressEvent,
                                this,
                                SLOT(DeselectInsertedVertex(vtkObject *)));
    this->ResetTrimBody();
  }
  this->Render();
}

void Viewer::ResetTrimBody() {
  if (trimmed_actor_)
    trimmed_actor_->GetProperty()->SetColor(snake_color_);
  trimmed_actor_ = NULL;
  trimmed_snake_ = NULL;
  is_trim_body_second_click_ = false;
  trim_body_index1_ = trim_body_index2_ = kBigNumber;
  renderer_->RemoveActor(on_snake_sphere1_);
  renderer_->RemoveActor(on_snake_sphere2_);
  renderer_->RemoveActor(off_snake_sphere_);
  inserted_point_.Fill(-1);
}

void Viewer::SelectBodyVertex() {
  picker_->Pick(qvtk_->GetInteractor()->GetEventPosition()[0],
                qvtk_->GetInteractor()->GetEventPosition()[1],
                0, renderer_);
  vtkActor *actor = picker_->GetActor();
  if (actor) {
    ActorSnakeMap::const_iterator it = actor_snake_map_.find(actor);
    if (it != actor_snake_map_.end()) {
      unsigned index = static_cast<unsigned>(picker_->GetPointId());

      if (is_trim_body_second_click_) {
        if (it->first == trimmed_actor_) {
          this->SetupSphere(it->second->GetPoint(index), on_snake_sphere2_,
                            sphere_color_);
          trim_body_index2_ = index;
          is_trim_body_second_click_ = false;
          trimmed_actor_->GetProperty()->SetColor(selected_snake_color_);
        }
      } else {
        this->SetupSphere(it->second->GetPoint(index), on_snake_sphere1_,
                          sphere_color_);
        if (trimmed_actor_) {
          if (it->first == trimmed_actor_) {
            trimmed_actor_->GetProperty()->SetColor(selected_snake_color_);
          } else {
            trimmed_actor_->GetProperty()->SetColor(snake_color_);
          }
        }
        trimmed_actor_ = it->first;
        trimmed_snake_ = it->second;
        trim_body_index1_ = index;
        is_trim_body_second_click_ = true;
      }
    }
  }
}

void Viewer::DeselectBodyVertex() {
  picker_->Pick(qvtk_->GetInteractor()->GetEventPosition()[0],
                qvtk_->GetInteractor()->GetEventPosition()[1],
                0, renderer_);
  vtkActor *actor = picker_->GetActor();
  if (actor) {
    if (actor == on_snake_sphere1_) {
      renderer_->RemoveActor(on_snake_sphere1_);
      trim_body_index1_ = kBigNumber;
      is_trim_body_second_click_ = false;
      trimmed_actor_->GetProperty()->SetColor(snake_color_);
    } else if (actor == on_snake_sphere2_) {
      renderer_->RemoveActor(on_snake_sphere2_);
      trim_body_index2_ = kBigNumber;
      is_trim_body_second_click_ = true;
      trimmed_actor_->GetProperty()->SetColor(snake_color_);
    }
  }
}

void Viewer::SelectInsertedVertex(vtkObject *obj) {
  vtkImagePlaneWidget *p = vtkImagePlaneWidget::SafeDownCast(obj);
  double point[3];
  p->GetCurrentCursorPosition(point);

  if (point[0] >= 0 && trim_body_index1_ != trim_body_index2_) {
    inserted_point_[0] = point[0];
    inserted_point_[1] = point[1];
    inserted_point_[2] = point[2];

    this->SetupSphere(inserted_point_, off_snake_sphere_, sphere_color_);
  }
}


void Viewer::DeselectInsertedVertex(vtkObject *obj) {
  picker_->Pick(qvtk_->GetInteractor()->GetEventPosition()[0],
                qvtk_->GetInteractor()->GetEventPosition()[1],
                0, renderer_);
  vtkActor *actor = picker_->GetActor();
  if (actor == off_snake_sphere_) {
    renderer_->RemoveActor(off_snake_sphere_);
    inserted_point_.Fill(-1);
  }
}

void Viewer::TrimBody() {
  if (!trimmed_snake_ || inserted_point_[0] < 0) return;
  trimmed_snake_->TrimAndInsert(trim_body_index1_, trim_body_index2_,
                                inserted_point_);
  trimmed_snake_->Resample();
  this->SetupSnake(trimmed_snake_, 0);
  inserted_point_.Fill(-1.0);
  trim_body_index1_ = trim_body_index2_ = kBigNumber;
  trimmed_actor_ = NULL;

  renderer_->RemoveActor(on_snake_sphere1_);
  renderer_->RemoveActor(on_snake_sphere2_);
  renderer_->RemoveActor(off_snake_sphere_);
}

void Viewer::ToggleDeleteJunction(bool state) {
  if (state) {
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::LeftButtonPressEvent,
                             this, SLOT(SelectJunction()));
    slot_connector_->Connect(qvtk_->GetInteractor(),
                             vtkCommand::RightButtonPressEvent,
                             this, SLOT(DeselectJunction()));
  } else {
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::LeftButtonPressEvent,
                                this, SLOT(SelectJunction()));
    slot_connector_->Disconnect(qvtk_->GetInteractor(),
                                vtkCommand::RightButtonPressEvent,
                                this, SLOT(DeselectJunction()));
  }
  this->Render();
}

void Viewer::SelectJunction() {
  picker_->Pick(qvtk_->GetInteractor()->GetEventPosition()[0],
                qvtk_->GetInteractor()->GetEventPosition()[1],
                0, renderer_);
  vtkActor *actor = picker_->GetActor();
  if (actor) {
    ActorPointMap::const_iterator it = actor_junctions_.find(actor);
    if (it != actor_junctions_.end()) {
      actor->GetProperty()->SetColor(selected_junction_color_);
      selected_junctions_.insert(*it);
    }
  }
}

void Viewer::DeselectJunction() {
  picker_->Pick(qvtk_->GetInteractor()->GetEventPosition()[0],
                qvtk_->GetInteractor()->GetEventPosition()[1],
                0, renderer_);
  vtkActor *actor = picker_->GetActor();
  if (actor) {
    ActorPointMap::iterator it = actor_junctions_.find(actor);
    if (it !=  actor_junctions_.end()) {
      actor->GetProperty()->SetColor(junction_color_);
      selected_junctions_.erase(actor);
    }
  }
}

void Viewer::RemoveSelectedJunctions(Junctions &junctions) {
  if (selected_junctions_.empty()) return;
  for (ActorPointMap::iterator it = selected_junctions_.begin();
       it != selected_junctions_.end(); ++it) {
    junctions.RemoveJunction(it->second);
    renderer_->RemoveActor(it->first);
    actor_junctions_.erase(it->first);
    it->first->Delete();
  }
  selected_junctions_.clear();
}

void Viewer::Render() {
  qvtk_->GetRenderWindow()->Render();
}

void Viewer::ResetCamera() {
  renderer_->ResetCamera();
}

void Viewer::LoadViewpoint(const std::string &filename) {
  std::ifstream infile(filename.c_str());
  std::string str;
  double x, y, z;

  infile >> str;
  if (str == "ClippingRange") {
    infile >> x >> y;
    camera_->SetClippingRange(x, y);
  }

  infile >> str;
  if (str == "CameraPosition") {
    infile >> x >> y >> z;
    camera_->SetPosition(x, y, z);
  }

  infile >> str;
  if (str == "CameraFocalPoint") {
    infile >> x >> y >> z;
    camera_->SetFocalPoint(x, y, z);
  }

  infile >> str;
  if (str == "CameraViewUp") {
    infile >> x >> y >> z;
    camera_->SetViewUp(x, y, z);
  }

  infile.close();
  this->Render();
}

void Viewer::SaveViewpoint(const std::string &filename) const {
  std::ofstream outfile(filename.c_str());
  double *pval;
  pval = camera_->GetClippingRange();
  outfile << "ClippingRange" << std::endl;
  outfile << pval[0] << " " << pval[1] << std::endl;

  pval = camera_->GetPosition();
  outfile << "CameraPosition" << std::endl;
  outfile << pval[0] << " " << pval[1] << " " << pval[2] << std::endl;

  pval = camera_->GetFocalPoint();
  outfile << "CameraFocalPoint" << std::endl;
  outfile << pval[0] << " " << pval[1] << " " << pval[2] << std::endl;

  pval = camera_->GetViewUp();
  outfile << "CameraViewUp" << std::endl;
  outfile << pval[0] << " " << pval[1] << " " << pval[2] << std::endl;

  outfile.close();
}


void Viewer::PrintScreenAsPNGImage(const std::string &filename) const {
  vtkNew<vtkWindowToImageFilter> filter;
  filter->SetInput(qvtk_->GetRenderWindow());
  // filter->SetMagnification(2);
  filter->SetInputBufferTypeToRGBA();
  filter->Update();
  vtkNew<vtkPNGWriter> png_writer;
  png_writer->SetInputConnection(filter->GetOutputPort());
  png_writer->SetFileName(filename.c_str());
  png_writer->Write();
}

void Viewer::PrintScreenAsTIFFImage(const std::string &filename) const {
  vtkNew<vtkWindowToImageFilter> filter;
  filter->SetInput(qvtk_->GetRenderWindow());
  // filter->SetMagnification(2);
  filter->SetInputBufferTypeToRGBA();
  filter->Update();

  std::string name = filename;
  name.replace(name.length()-3, 3, "tif");

  vtkNew<vtkTIFFWriter> tiff_writer;
  tiff_writer->SetInputConnection(filter->GetOutputPort());
  tiff_writer->SetFileName(name.c_str());
  tiff_writer->Write();
}

void Viewer::PrintScreenAsVectorImage(const std::string &filename) const {
  std::string::size_type dot_position = filename.find_last_of(".");
  std::string name = filename;
  if (dot_position == std::string::npos) {
    std::cerr << "No suffix specified! Abort." << std::endl;
    return;
  } else {
    name = filename.substr(0, dot_position);
  }

  vtkNew<vtkGL2PSExporter> exp;
  exp->SetRenderWindow(qvtk_->GetRenderWindow());
  exp->CompressOff();
  exp->SetSortToBSP();
  exp->SetFilePrefix(name.c_str());
  exp->SetTitle("SOAX Screenshot");
  exp->TextAsPathOn();
  exp->SetFileFormatToEPS();
  exp->Write();
}

}  // namespace soax
