#include "viewer.h"
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

#include "itkStatisticsImageFilter.h"


namespace soax {

const double Viewer::kWhite[3] = {1.0, 1.0, 1.0};
const double Viewer::kGray[3] = {0.8, 0.8, 0.8};
const double Viewer::kRed[3] = {1.0, 0.0, 0.0};
const double Viewer::kMagenta[3] = {1.0, 0.0, 1.0};
const double Viewer::kYellow[3] = {1.0, 1.0, 0.0};
const double Viewer::kGreen[3] = {0.0, 1.0, 0.0};
const double Viewer::kCyan[3] = {0.0, 1.0, 1.0};
const double Viewer::kBlue[3] = {0.0, 0.0, 1.0};

Viewer::Viewer() {
  qvtk_ = new QVTKWidget;
  renderer_ = vtkSmartPointer<vtkRenderer>::New();
  qvtk_->GetRenderWindow()->AddRenderer(renderer_);
  camera_ = vtkSmartPointer<vtkCamera>::New();
  renderer_->SetActiveCamera(camera_);

  for (unsigned i = 0; i < kDimension; i++) {
    slice_planes_[i] = vtkSmartPointer<vtkImagePlaneWidget>::New();
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
}

Viewer::~Viewer() {
  delete qvtk_;
}

void Viewer::SetupImage(ImageType::Pointer image) {
  typedef itk::ImageToVTKImageFilter<ImageType>  ConnectorType;
  ConnectorType::Pointer connector = ConnectorType::New();
  connector->SetInput(image);
  connector->Update();
  // vtkSmartPointer<vtkImageCast> caster = vtkSmartPointer<vtkImageCast>::New();
  // caster->SetInputData(connector->GetOutput());
  // caster->SetOutputScalarTypeToUnsignedShort();
  // caster->Update();

  typedef itk::StatisticsImageFilter<ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  filter->Update();

  this->SetupSlicePlanes(connector->GetOutput(), filter->GetMinimum(),
                         filter->GetMaximum());
  // this->SetupMIPRendering(caster->GetOutput());
}

void Viewer::SetupSlicePlanes(vtkImageData *data, double min_intensity,
                              double max_intensity) {
  double window = max_intensity - min_intensity;
  double level = min_intensity + window/2;

  for (unsigned i = 0; i < kDimension; ++i) {
    slice_planes_[i]->SetInputData(data);
    slice_planes_[i]->SetPlaneOrientation(i);
    slice_planes_[i]->UpdatePlacement();
    slice_planes_[i]->SetWindowLevel(window, level);
  }
}

void Viewer::ToggleSlicePlanes(bool state) {
  if (state)
    this->ShowSlicePlanes();
  else
    this->HideSlicePlanes();
  this->Render();
}

void Viewer::ShowSlicePlanes() {
  for (unsigned i = 0; i < kDimension; ++i) {
    slice_planes_[i]->On();
  }
  renderer_->ResetCamera();
}

void Viewer::HideSlicePlanes() {
  for (unsigned i = 0; i < kDimension; ++i) {
    slice_planes_[i]->Off();
  }
}

void Viewer::Render() {
  qvtk_->GetRenderWindow()->Render();
}

} // namespace soax
