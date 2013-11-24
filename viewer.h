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
class vtkPolyData;


namespace soax {

class Viewer : public QObject {
  Q_OBJECT

 public:
  Viewer();
  ~Viewer();
  void Reset();
  QVTKWidget *qvtk() const {return qvtk_;}

  void SetupImage(ImageType::Pointer image);
  void SetupSnakes(const SnakeContainer &snakes, unsigned category = 0);
  void RemoveSnakes();
  void SetupJunctions(const PointContainer &points);
  void RemoveJunctions();
  void SetupUpperRightCornerText();
  void set_snake_filename(const std::string &filename) {
    snake_filename_ = filename;
  }
  void set_comapring_snake_filename1(const std::string &filename) {
    comparing_snake_filename1_ = filename;
  }
  void set_comapring_snake_filename2(const std::string &filename) {
    comparing_snake_filename2_ = filename;
  }

  void Render();
  void ResetCamera();

  void LoadViewpoint(const std::string &filename);
  void SaveViewpoint(const std::string &filename) const;
  void PrintScreenAsPNGImage(const std::string &filename) const;
  void PrintScreenAsTIFFImage(const std::string &filename) const;
  void PrintScreenAsVectorImage(const std::string &filename) const;

 public slots:
  void ToggleSlicePlanes(bool state);
  void ToggleMIPRendering(bool state);
  void ToggleOrientationMarker(bool state);
  void ToggleCornerText(bool state);
  void ToggleBoundingBox(bool state);
  void ToggleCubeAxes(bool state);
  void ToggleSnakes(bool state);
  void ToggleJunctions(bool state);

 private:
  typedef std::map<vtkActor *, Snake *> ActorSnakeMap;
  typedef std::map<Snake *, vtkActor *> SnakeActorMap;
  typedef std::map<vtkActor *, PointType> ActorPointMap;


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
  void SetupSnake(Snake *snake, unsigned category);
  vtkActor * ActSnake(Snake *snake);
  vtkActor * ActSnakeCell(Snake *snake, unsigned start, unsigned end);
  vtkPolyData * MakePolyData(Snake *snake,unsigned start, unsigned end);
  void SetupEvolvingActorProperty(vtkActor *actor);
  void SetupComparingActorProperty(vtkActor *actor);
  void SetupAnotherComparingActorProperty(vtkActor *actor);
  void SetupSphere(const PointType &point, vtkActor *sphere);


  QVTKWidget *qvtk_;
  vtkSmartPointer<vtkRenderer> renderer_;
  vtkSmartPointer<vtkCamera> camera_;
  vtkSmartPointer<vtkImagePlaneWidget> slice_planes_[kDimension];
  vtkSmartPointer<vtkVolume> volume_;
  vtkSmartPointer<vtkOrientationMarkerWidget> orientation_marker_;
  vtkSmartPointer<vtkCornerAnnotation> corner_text_;
  vtkSmartPointer<vtkCubeAxesActor> cube_axes_;
  vtkSmartPointer<vtkActor> bounding_box_;

  ActorSnakeMap actor_snakes_;
  SnakeActorMap snake_actors_;
  ActorSnakeMap actor_selected_snakes_;
  ActorPointMap junctions_;
  ActorPointMap selected_junctions_;

  double snake_opacity_;
  double comparing_snakes1_opacity_;
  double comparing_snakes2_opacity_;

  double *snake_color_;
  double *comparing_snakes1_color_;
  double *comparing_snakes2_color_;

  double snake_width_;
  double comparing_snakes1_width_;
  double comparing_snakes2_width_;

  double junction_radius_;
  double *junction_color_;

  std::string snake_filename_;
  std::string comparing_snake_filename1_;
  std::string comparing_snake_filename2_;


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
