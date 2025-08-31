/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file defines the visulization class for SOAX.
 */


#ifndef VIEWER_H_
#define VIEWER_H_

#include <string>
#include <vector>
#include <map>
#include <QObject>  // NOLINT(build/include_order)
#include "vtkSmartPointer.h"
#include "./global.h"
#include "./junctions.h"
#include <QVTKOpenGLWidget.h>
#include "vtkEventQtSlotConnect.h"

class QVTKOpenGLWidget;
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
class vtkObject;
class vtkPointPicker;

namespace soax {

class Viewer : public QObject {
  Q_OBJECT

public:
  Viewer(QVTKOpenGLWidget *qvtk = nullptr);
  ~Viewer();
  void Reset();
  QVTKOpenGLWidget *qvtk() const {return qvtk_;}

  void SetupImage(ImageType::Pointer image);

  double window() const {return window_;}
  void set_window(double window) {window_ = window;}

  double level() const {return level_;}
  void set_level(double level) {level_ = level;}

  double mip_min_intensity() const {return mip_min_intensity_;}
  void set_mip_min_intensity(double mip_min) {
    mip_min_intensity_ = mip_min;
  }

  double mip_max_intensity() const {return mip_max_intensity_;}
  void set_mip_max_intensity(double mip_max) {
    mip_max_intensity_ = mip_max;
  }

  double clip_span() const {return clip_span_;}
  void set_clip_span(double span) {clip_span_ = span;}

  unsigned color_segment_step() const {return color_segment_step_;}
  void set_color_segment_step(unsigned step) {color_segment_step_ = step;}

  Snake * trimmed_snake() const {return trimmed_snake_;}
  Snake * selected_snake() const {return selected_snake_;}

  void UpdateWindowLevel(double window, double level);
  void UpdateMIPIntensityRange(double min, double max);

  void SetupSnakesAsOneActor(const SnakeContainer &snakes);
  void SetupSnakes(const SnakeContainer &snakes, unsigned category = 0);
  void SetupSnake(Snake *snake, unsigned category);
  void ChangeSnakeColor(Snake *s, double *color);
  void RemoveSnakes();
  void RemoveSnake(Snake *s);
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

  const SnakeSet &selected_snakes() const {return selected_snakes_;}
  void RemoveSelectedSnakes();
  void TrimTip();
  void ExtendTip();
  void TrimBody();
  void RemoveSelectedJunctions(Junctions &junctions);

  static double *Gray() {return kGray;}
  static double *Red() {return kRed;}
  static double *Magenta() {return kMagenta;}
  static double *Yellow() {return kYellow;}
  static double *Green() {return kGreen;}
  static double *Cyan() {return kCyan;}
  static double *Blue() {return kBlue;}

 public slots:  // NOLINT(whitespace/indent)
  void ToggleSlicePlanes(bool state);
  void ToggleMIPRendering(bool state);
  void ToggleOrientationMarker(bool state);
  void ToggleCornerText(bool state);
  void ToggleBoundingBox(bool state);
  void ToggleCubeAxes(bool state);
  void ToggleSnakes(bool state);
  void ToggleJunctions(bool state);
  void ToggleClipSnakes(bool state);
  void ColorByAzimuthalAngle(bool state);
  void ColorByPolarAngle(bool state);

  void ToggleNone(bool state);
  void ToggleDeleteSnake(bool state);
  void ToggleTrimTip(bool state);
  void ToggleExtendTip(bool state);
  void ToggleTrimBody(bool state);
  void ToggleDeleteJunction(bool state);

 private slots:  // NOLINT(whitespace/indent)
  void SetupClippedSnakes(vtkObject *obj);
  void SelectSnakeForView();
  void DeselectSnakeForView();
  void SelectSnakeForDeletion();
  void DeselectSnakeForDeletion();
  void SelectVertex();
  void DeselectVertex();
  void SelectExtendVertex(vtkObject *obj);
  void DeselectExtendVertex(vtkObject *obj);
  void SelectBodyVertex();
  void DeselectBodyVertex();
  void SelectInsertedVertex(vtkObject *obj);
  void DeselectInsertedVertex(vtkObject *obj);
  void SelectJunction();
  void DeselectJunction();

 private:
  typedef std::map<vtkActor *, Snake *> ActorSnakeMap;
  typedef std::map<Snake *, vtkActor *> SnakeActorMap;
  typedef std::map<vtkActor *, PointType> ActorPointMap;
  typedef std::vector<vtkActor *> ActorContainer;

  void SetupSlicePlanes(vtkSmartPointer<vtkImageData> data);
  void SetupMIPRendering(vtkSmartPointer<vtkImageData> data);
  void SetupOrientationMarker();
  void SetupUpperLeftCornerText(unsigned min_intensity,
                                unsigned max_intensity);
  void SetupBoundingBox(vtkSmartPointer<vtkVolume> volume);
  void SetupCubeAxes(vtkSmartPointer<vtkImageData> image);

  void UpdateJunctionRadius(ImageType::Pointer image);

  vtkActor * ActSnake(Snake *snake);
  vtkActor * ActSnakeSegments(Snake *snake, unsigned start, unsigned end);
  vtkPolyData * MakePolyDataForMultipleSnakes(const SnakeContainer &snakes);
  vtkPolyData * MakePolyData(Snake *snake, unsigned start, unsigned end);
  void SetupEvolvingActorProperty(vtkActor *actor);
  void SetupComparingActorProperty(vtkActor *actor);
  void SetupAnotherComparingActorProperty(vtkActor *actor);
  void SetupSphere(const PointType &point, vtkActor *sphere, double *color);

  vtkPolyData * MakeClippedPolyData(unsigned axis, double position);

  void ColorSnakes(bool state, bool azimuthal);
  void SetupColorSegments(bool azimuthal);
  void ComputeColor(const VectorType &vector, bool azimuthal, double *color);
  void ComputeThetaPhi(const VectorType &vector, double *theta, double *phi);
  void ComputeRGBFromHue(double hue, double *color);
  void RemoveColorSegments();
  void ResetTrimTip();
  void ResetExtendTip();
  void ResetTrimBody();


  QVTKOpenGLWidget *qvtk_;
  vtkGenericOpenGLRenderWindow *render_window_;  
  vtkSmartPointer<vtkRenderer> renderer_;
  vtkSmartPointer<vtkCamera> camera_;
  vtkImagePlaneWidget *slice_planes_[kDimension];
  vtkSmartPointer<vtkVolume> volume_;
  vtkSmartPointer<vtkOrientationMarkerWidget> orientation_marker_;
  vtkSmartPointer<vtkCornerAnnotation> corner_text_;
  vtkSmartPointer<vtkCubeAxesActor> cube_axes_;
  vtkSmartPointer<vtkActor> bounding_box_;
  vtkSmartPointer<vtkEventQtSlotConnect> slot_connector_;
  vtkSmartPointer<vtkActor> clipped_actor_;
  ActorContainer color_segments_;
  vtkSmartPointer<vtkPointPicker> picker_;

  double window_;
  double level_;
  double mip_min_intensity_;
  double mip_max_intensity_;
  double clip_span_;
  unsigned color_segment_step_;

  ActorSnakeMap actor_snake_map_;
  SnakeActorMap snake_actor_map_;
  ActorPointMap actor_junctions_;
  SnakeSet selected_snakes_;
  ActorPointMap selected_junctions_;

  vtkActor *snakes_actor_;
  vtkActor *on_snake_sphere1_;
  vtkActor *on_snake_sphere2_;
  vtkActor *off_snake_sphere_;
  vtkActor *trimmed_actor_;
  Snake *selected_snake_;
  Snake *trimmed_snake_;
  unsigned trim_tip_index_;
  unsigned trim_body_index1_;
  unsigned trim_body_index2_;
  PointType inserted_point_;
  bool is_trim_body_second_click_;

  double snake_opacity_;
  double comparing_snakes1_opacity_;
  double comparing_snakes2_opacity_;

  double *snake_color_;
  double *comparing_snakes1_color_;
  double *comparing_snakes2_color_;
  double *selected_snake_color_;
  double *sphere_color_;

  double snake_width_;
  double comparing_snakes1_width_;
  double comparing_snakes2_width_;

  double junction_radius_;
  double *junction_color_;
  double *selected_junction_color_;

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

}  // namespace soax

#endif  // VIEWER_H_
