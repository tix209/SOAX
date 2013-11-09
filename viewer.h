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
  void ToggleCornerText(bool state);
  void ToggleBoundingBox(bool state);
  void ToggleCubeAxes(bool state);

 private:
  void SetupSlicePlanes(vtkImageData *data, double min_intensity,
                        double max_intensity);
  void SetupMIPRendering(vtkImageData *data,
                                 double min_intensity,
                                 double max_intensity);
  void SetupOrientationMarker();
  void SetupUpperLeftCornerText(unsigned min_intensity,
                                unsigned max_intensity);
  void SetupBoundingBox();
  void SetupCubeAxes(vtkImageData *image);

  QVTKWidget *qvtk_;
  vtkSmartPointer<vtkRenderer> renderer_;
  vtkSmartPointer<vtkCamera> camera_;
  vtkSmartPointer<vtkImagePlaneWidget> slice_planes_[kDimension];
  vtkSmartPointer<vtkVolume> volume_;
  vtkSmartPointer<vtkOrientationMarkerWidget> orientation_marker_;
  vtkSmartPointer<vtkCornerAnnotation> corner_text_;
  vtkSmartPointer<vtkCubeAxesActor> cube_axes_;
  vtkSmartPointer<vtkActor> bounding_box_;


  static double kWhite[3];
  static double kGray[3];
  static double kRed[3];
  static double kMagenta[3];
  static double kYellow[3];
  static double kGreen[3];
  static double kCyan[3];
  static double kBlue[3];

  DISALLOW_COPY_AND_ASSIGN(Viewer);
};

} // namespace soax

#endif // SOAX_VIEWER_H_
