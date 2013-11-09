#ifndef SOAX_VIEWER_H_
#define SOAX_VIEWER_H_

#include <QObject>
#include "vtkSmartPointer.h"
#include "global.h"

class QVTKWidget;
class vtkImagePlaneWidget;
class vtkRenderer;
class vtkCamera;
class vtkVolume;
class vtkOrientationMarkerWidget;
class vtkCornerAnnotation;
class vtkCubeAxesActor;
class vtkActor;
class vtkImageData;

namespace soax {

class Viewer : public QObject {
  Q_OBJECT

 public:
  Viewer();
  ~Viewer();

  QVTKWidget *qvtk() const {return qvtk_;}

  void SetupImage(ImageType::Pointer image);

  void Render();

 public slots:
  void ToggleSlicePlanes(bool state);
  void ToggleMIPRendering(bool state);
  void ToggleOrientationMarker(bool state);
  void ToggleScreenInformation(bool state);
  void ToggleBoundingBox(bool state);
  void ToggleCubeAxes(bool state);

 private:
  void SetupSlicePlanes(vtkImageData *data, double min_intensity,
                        double max_intensity);
  void SetupMIPRendering(vtkImageData *data,
                                 double min_intensity,
                                 double max_intensity);
  void SetupOrientationMarker();

  QVTKWidget *qvtk_;
  vtkSmartPointer<vtkRenderer> renderer_;
  vtkSmartPointer<vtkCamera> camera_;
  vtkSmartPointer<vtkImagePlaneWidget> slice_planes_[kDimension];
  vtkSmartPointer<vtkVolume> volume_;
  vtkSmartPointer<vtkOrientationMarkerWidget> orientation_marker_;
  vtkSmartPointer<vtkCornerAnnotation> corner_text_;
  vtkSmartPointer<vtkCubeAxesActor> cube_axes_;
  vtkSmartPointer<vtkActor> bounding_box_;


  static const double kWhite[3];
  static const double kGray[3];
  static const double kRed[3];
  static const double kMagenta[3];
  static const double kYellow[3];
  static const double kGreen[3];
  static const double kCyan[3];
  static const double kBlue[3];

  DISALLOW_COPY_AND_ASSIGN(Viewer);
};

} // namespace soax

#endif // SOAX_VIEWER_H_
