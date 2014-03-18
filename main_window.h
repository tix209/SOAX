#ifndef SOAX_MAINWINDOW_H_
#define SOAX_MAINWINDOW_H_

#include <QMainWindow>
#include "global.h"

class QProgressBar;
class QActionGroup;


namespace soax {

class Multisnake;
class Viewer;
class ParametersDialog;
class AnalysisOptionsDialog;
class ViewOptionsDialog;

/*
 * SOAX GUI main window.
 */
class MainWindow : public QMainWindow {
  Q_OBJECT

 public:
  MainWindow();
  ~MainWindow();

  // Disable right click on the toolbar
  virtual QMenu *createPopupMenu() {return NULL;}

 private slots:
  void OpenImage();
  void SaveAsIsotropicImage();
  void LoadParameters();
  void SaveParameters();
  void LoadSnakes();
  void SaveSnakes();
  void CompareSnakes();
  void CompareAnotherSnakes();
  void LoadJFilamentSnakes();
  void SaveJFilamentSnakes();
  void CloseSession();

  void EditSnake();

  // void ToggleSnakeDisplay(bool state);
  void ShowViewOptions();

  void InitializeSnakes();
  void DeformSnakes();
  void DeformSnakesInAction();
  void DeformOneSnake();
  void CutSnakes();
  void GroupSnakes();

  void ComputeSphericalOrientation();
  void ComputeRadialOrientation();
  void ComputePointDensity();
  void ComputeCurvature();
  void ShowAnalysisOptions();

  void ShowParametersDialog();
  void LoadViewpoint();
  void SaveViewpoint();
  void SaveSnapshot();

  void AboutSOAX();

 private:
  void CreateActions();
  void CreateFileMenuActions();
  void CreateEditMenuActions();
  void CreateViewMenuActions();
  void CreateProcessMenuActions();
  void CreateAnalysisMenuActions();
  void CreateToolsMenuActions();
  void CreateHelpMenuActions();

  void CreateMenus();
  void CreateToolBar();
  void CreateStatusBar();

  void ResetActions();
  QString GetLastDirectory(const std::string &filename) const;
  void SetParameters();

  QWidget *central_widget_;
  ParametersDialog *parameters_dialog_;
  ViewOptionsDialog *view_options_dialog_;
  AnalysisOptionsDialog *analysis_options_dialog_;

  // Actions in File menu
  QAction *open_image_;
  QAction *save_as_isotropic_image_;
  QAction *load_parameters_;
  QAction *save_parameters_;
  QAction *load_snakes_;
  QAction *save_snakes_;
  QAction *load_jfilament_snakes_;
  QAction *save_jfilament_snakes_;
  QAction *compare_snakes_;
  QAction *compare_another_snakes_;
  QAction *close_session_;
  QAction *exit_;

  // Actions in Edit menu
  QAction *toggle_none_;
  QAction *toggle_delete_snake_;
  QAction *toggle_trim_tip_;
  QAction *toggle_extend_tip_;
  QAction *toggle_trim_body_;
  QAction *toggle_delete_junction_;
  QAction *edit_snake_;
  QActionGroup *snake_edit_group_;

  // Actions in View menu
  QAction *toggle_planes_;
  QAction *toggle_mip_;
  QAction *toggle_orientation_marker_;
  QAction *toggle_corner_text_;
  QAction *toggle_bounding_box_;
  QAction *toggle_cube_axes_;
  QAction *toggle_snakes_;
  QAction *toggle_junctions_;
  QAction *toggle_clip_;
  QAction *toggle_color_azimuthal_;
  QAction *toggle_color_polar_;
  QAction *show_view_options_;
  QActionGroup *snake_view_group_;

  // Actions in Process menu
  QAction *initialize_snakes_;
  QAction *deform_snakes_;
  QAction *deform_snakes_in_action_;
  QAction *deform_one_snake_;
  QAction *cut_snakes_;
  QAction *group_snakes_;

  // Actions in Analysis menu
  QAction *compute_spherical_orientation_;
  QAction *compute_radial_orientation_;
  QAction *compute_point_density_;
  QAction *compute_curvature_;
  QAction *show_analysis_options_;

  // Actions in Tools menu
  QAction *show_parameters_;
  QAction *load_viewpoint_;
  QAction *save_viewpoint_;
  QAction *save_snapshot_;

  // Actions in Help menu
  QAction *about_soax_;
  QAction *about_qt_;

  // Menus
  QMenu *file_;
  QMenu *edit_;
  QMenu *snake_edit_submenu_;
  QMenu *view_;
  QMenu *process_;
  QMenu *analysis_;
  QMenu *actin_cable_submenu_;
  QMenu *contractile_ring_submenu_;
  QMenu *fibrin_submenu_;
  QMenu *droplet_submenu_;
  QMenu *tools_;
  QMenu *help_;

  QToolBar *toolbar_;
  QProgressBar *progress_bar_;
  int message_timeout_; // in milliseconds

  // Complete image file path and name
  std::string image_filename_;
  std::string parameter_filename_;
  std::string snake_filename_;
  std::string viewpoint_filename_;
  std::string snapshot_filename_;

  Multisnake *multisnake_;
  Viewer *viewer_;

  DISALLOW_COPY_AND_ASSIGN(MainWindow);
};

} // namespace soax

#endif // SOAX_MAINWINDOW_H_
