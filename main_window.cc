/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file implements the main window class for SOAX.
 */

#include <QtWidgets>
#include "./main_window.h"
#include "./multisnake.h"
#include "./viewer.h"
#include "./parameters_dialog.h"
#include "./view_options_dialog.h"
#include "./solver_bank.h"
#include "./analysis_options_dialog.h"
#include "./utility.h"


namespace soax {

MainWindow::MainWindow() : message_timeout_(0) {
  central_widget_ = new QWidget(this);
  multisnake_ = new Multisnake;
  QVTKOpenGLWidget* qvtk_widget = new QVTKOpenGLWidget(this);
  viewer_ = new Viewer(qvtk_widget);
  parameters_dialog_ = new ParametersDialog(this);
  view_options_dialog_ = new ViewOptionsDialog(this);
  analysis_options_dialog_ = new AnalysisOptionsDialog(this);

  QVBoxLayout* layout = new QVBoxLayout;
  layout->addWidget(viewer_->qvtk());
  central_widget_->setLayout(layout);
  setCentralWidget(central_widget_);
  setWindowIcon(QIcon(":/icon/letter-x.png"));
  setWindowTitle("SOAX");
  setUnifiedTitleAndToolBarOnMac(true);

  this->CreateActions();
  this->CreateMenus();
  this->CreateToolBar();
  this->CreateStatusBar();
  this->ResetActions();
}


MainWindow::~MainWindow() {
  delete multisnake_;
  delete viewer_;
}

void MainWindow::CreateActions() {
  this->CreateFileMenuActions();
  this->CreateEditMenuActions();
  this->CreateViewMenuActions();
  this->CreateProcessMenuActions();
  this->CreateAnalysisMenuActions();
  this->CreateToolsMenuActions();
  this->CreateHelpMenuActions();
}

void MainWindow::CreateFileMenuActions() {
  open_image_ = new QAction(tr("&Open Image"), this);
  open_image_->setShortcut(QKeySequence::Open);
  open_image_->setIcon(QIcon(":/icon/Open.png"));
  connect(open_image_, SIGNAL(triggered()),
          this, SLOT(OpenImage()));

  save_as_isotropic_image_ = new QAction(tr("Save as Isotropic Image"),
                                         this);
  connect(save_as_isotropic_image_, SIGNAL(triggered()),
          this, SLOT(SaveAsIsotropicImage()));

  load_parameters_ = new QAction(tr("Load Pa&rameters"), this);
  load_parameters_->setShortcut(Qt::CTRL + Qt::Key_R);
  load_parameters_->setIcon(QIcon(":/icon/Properties.png"));
  connect(load_parameters_, SIGNAL(triggered()),
          this, SLOT(LoadParameters()));

  save_parameters_ = new QAction(tr("Save Parameters"), this);
  connect(save_parameters_, SIGNAL(triggered()),
          this, SLOT(SaveParameters()));

  load_snakes_ = new QAction(tr("Load Snakes"), this);
  load_snakes_->setShortcut(Qt::CTRL + Qt::Key_L);
  connect(load_snakes_, SIGNAL(triggered()), this, SLOT(LoadSnakes()));

  save_snakes_ = new QAction(tr("Save Snakes"), this);
  save_snakes_->setShortcut(QKeySequence::Save);
  save_snakes_->setIcon(QIcon(":/icon/Save.png"));
  connect(save_snakes_, SIGNAL(triggered()), this, SLOT(SaveSnakes()));

  load_jfilament_snakes_ = new QAction(tr("Load JFilament Snakes"), this);
  connect(load_jfilament_snakes_, SIGNAL(triggered()),
          this, SLOT(LoadJFilamentSnakes()));

  save_jfilament_snakes_ = new QAction(tr("Save JFilament Snakes"), this);
  connect(save_jfilament_snakes_, SIGNAL(triggered()),
          this, SLOT(SaveJFilamentSnakes()));

  compare_snakes_ = new QAction(tr("Compare Snakes"), this);
  compare_snakes_->setShortcut(Qt::CTRL + Qt::Key_C);
  connect(compare_snakes_, SIGNAL(triggered()),
          this, SLOT(CompareSnakes()));

  compare_another_snakes_ = new QAction(tr("Compare Another Snakes"), this);
  connect(compare_another_snakes_, SIGNAL(triggered()),
          this, SLOT(CompareAnotherSnakes()));

  close_session_ = new QAction(tr("Close Session"), this);
  close_session_->setIcon(QIcon(":/icon/Logout.png"));
  close_session_->setShortcut(QKeySequence::Close);
  connect(close_session_, SIGNAL(triggered()), this, SLOT(CloseSession()));

  exit_ = new QAction(tr("E&xit"), this);
  exit_->setShortcut(QKeySequence::Quit);
  connect(exit_, SIGNAL(triggered()), this, SLOT(close()));
}

void MainWindow::CreateEditMenuActions() {
  toggle_none_ = new QAction(tr("Normal Mode"), this);
  toggle_none_->setIcon(QIcon(":/icon/Cancel.png"));
  toggle_none_->setCheckable(true);
  connect(toggle_none_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleNone(bool)));

  toggle_delete_snake_ = new QAction(tr("Delete Snake Mode"), this);
  toggle_delete_snake_->setIcon(QIcon(":/icon/Delete.png"));
  toggle_delete_snake_->setCheckable(true);
  connect(toggle_delete_snake_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleDeleteSnake(bool)));

  toggle_trim_tip_ = new QAction(tr("Trim Tip Mode"), this);
  toggle_trim_tip_->setCheckable(true);
  toggle_trim_tip_->setIcon(QIcon(":/icon/letter-t.png"));
  connect(toggle_trim_tip_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleTrimTip(bool)));

  toggle_extend_tip_ = new QAction(tr("Extend Tip Mode"), this);
  toggle_extend_tip_->setCheckable(true);
  toggle_extend_tip_->setIcon(QIcon(":/icon/letter-e.png"));
  connect(toggle_extend_tip_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleExtendTip(bool)));

  toggle_trim_body_ = new QAction(tr("Trim Body Mode"), this);
  toggle_trim_body_->setCheckable(true);
  toggle_trim_body_->setIcon(QIcon(":/icon/letter-b.png"));
  connect(toggle_trim_body_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleTrimBody(bool)));

  toggle_delete_junction_ = new QAction(tr("Delete Junction Mode"), this);
  toggle_delete_junction_->setCheckable(true);
  toggle_delete_junction_->setIcon(QIcon(":/icon/Remove.png"));
  connect(toggle_delete_junction_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleDeleteJunction(bool)));

  edit_snake_ = new QAction(tr("Edit Snake"), this);
  edit_snake_->setShortcut(Qt::Key_Space);
  connect(edit_snake_, SIGNAL(triggered()), this, SLOT(EditSnake()));

  snake_edit_group_ = new QActionGroup(this);
  snake_edit_group_->addAction(toggle_none_);
  snake_edit_group_->addAction(toggle_delete_snake_);
  snake_edit_group_->addAction(toggle_trim_tip_);
  snake_edit_group_->addAction(toggle_extend_tip_);
  snake_edit_group_->addAction(toggle_trim_body_);
  snake_edit_group_->addAction(toggle_delete_junction_);
  snake_edit_group_->setExclusive(true);
  toggle_none_->setChecked(true);
}

void MainWindow::CreateViewMenuActions() {
  toggle_planes_ = new QAction(tr("Slice Planes"), this);
  toggle_planes_->setIcon(QIcon(":/icon/Picture.png"));
  toggle_planes_->setCheckable(true);
  connect(toggle_planes_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleSlicePlanes(bool)));

  toggle_mip_ = new QAction(tr("MIP Rendering"), this);
  toggle_mip_->setIcon(QIcon(":/icon/Globe.png"));
  toggle_mip_->setCheckable(true);
  connect(toggle_mip_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleMIPRendering(bool)));

  toggle_orientation_marker_ = new QAction(tr("Orientation Marker"), this);
  toggle_orientation_marker_->setCheckable(true);
  connect(toggle_orientation_marker_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleOrientationMarker(bool)));

  toggle_corner_text_ = new QAction(tr("Corner Texts"), this);
  toggle_corner_text_->setCheckable(true);
  connect(toggle_corner_text_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleCornerText(bool)));

  toggle_bounding_box_ = new QAction(tr("Bounding Box"), this);
  toggle_bounding_box_->setCheckable(true);
  connect(toggle_bounding_box_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleBoundingBox(bool)));

  toggle_cube_axes_ = new QAction(tr("Cube Axes"), this);
  toggle_cube_axes_->setCheckable(true);
  connect(toggle_cube_axes_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleCubeAxes(bool)));

  toggle_snakes_ = new QAction(tr("Snakes"), this);
  toggle_snakes_->setIcon(QIcon(":/icon/Synchronize.png"));
  toggle_snakes_->setCheckable(true);
  connect(toggle_snakes_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleSnakes(bool)));

  toggle_junctions_ = new QAction(tr("Junctions"), this);
  toggle_junctions_->setIcon(QIcon(":/icon/Positive.png"));
  toggle_junctions_->setCheckable(true);
  connect(toggle_junctions_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleJunctions(bool)));

  toggle_clip_ = new QAction(tr("Show Snakes Locally"), this);
  toggle_clip_->setCheckable(true);
  toggle_clip_->setIcon(QIcon(":/icon/Search.png"));
  connect(toggle_clip_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleClipSnakes(bool)));

  toggle_color_azimuthal_ = new QAction(
      tr("Color Snakes by Azimuthal Angle"), this);
  toggle_color_azimuthal_->setCheckable(true);
  connect(toggle_color_azimuthal_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ColorByAzimuthalAngle(bool)));

  toggle_color_polar_ = new QAction(tr("Color Snakes by Polar Angle"), this);
  toggle_color_polar_->setCheckable(true);
  connect(toggle_color_polar_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ColorByPolarAngle(bool)));

  show_view_options_ = new QAction(tr("Options"), this);
  connect(show_view_options_, SIGNAL(triggered()),
          this, SLOT(ShowViewOptions()));

  snake_view_group_ = new QActionGroup(this);
  snake_view_group_->addAction(toggle_snakes_);
  snake_view_group_->addAction(toggle_clip_);
  snake_view_group_->addAction(toggle_color_azimuthal_);
  snake_view_group_->addAction(toggle_color_polar_);
  snake_view_group_->setExclusive(true);
}

void MainWindow::CreateProcessMenuActions() {
  initialize_snakes_ = new QAction(tr("Initialize Snakes"), this);
  initialize_snakes_->setIcon(QIcon(":/icon/Add.png"));
  initialize_snakes_->setShortcut(Qt::CTRL + Qt::Key_A);
  connect(initialize_snakes_, SIGNAL(triggered()),
          this, SLOT(InitializeSnakes()));

  deform_snakes_ = new QAction(tr("Deform Snakes"), this);
  deform_snakes_->setIcon(QIcon(":/icon/Play.png"));
  deform_snakes_->setShortcut(Qt::CTRL + Qt::Key_D);
  connect(deform_snakes_, SIGNAL(triggered()), this, SLOT(DeformSnakes()));

  deform_snakes_in_action_ = new QAction(
      tr("Deform Snakes in Action"), this);
  deform_snakes_in_action_->setIcon(QIcon(":/icon/Refresh.png"));
  connect(deform_snakes_in_action_, SIGNAL(triggered()),
          this, SLOT(DeformSnakesInAction()));

  deform_one_snake_ = new QAction(tr("Deform One Snake"), this);
  deform_one_snake_->setIcon(QIcon(":/icon/Pause.png"));
  connect(deform_one_snake_, SIGNAL(triggered()),
          this, SLOT(DeformOneSnake()));

  cut_snakes_ = new QAction(tr("Cut Snakes at Junctions"), this);
  cut_snakes_->setIcon(QIcon(":/icon/Cut.png"));
  connect(cut_snakes_, SIGNAL(triggered()), this, SLOT(CutSnakes()));

  group_snakes_ = new QAction(tr("Group Snakes"), this);
  group_snakes_->setIcon(QIcon(":/icon/Favorites.png"));
  connect(group_snakes_, SIGNAL(triggered()), this, SLOT(GroupSnakes()));
}

void MainWindow::CreateAnalysisMenuActions() {
  compute_spherical_orientation_ = new QAction(
      tr("Compute Spherical Orientation"), this);
  connect(compute_spherical_orientation_, SIGNAL(triggered()),
          this, SLOT(ComputeSphericalOrientation()));

  compute_radial_orientation_ = new QAction(
      tr("Compute Radial Orientation"), this);
  connect(compute_radial_orientation_, SIGNAL(triggered()),
          this, SLOT(ComputeRadialOrientation()));

  compute_point_density_ = new QAction(tr("Compute Point Density"), this);
  connect(compute_point_density_, SIGNAL(triggered()),
          this, SLOT(ComputePointDensity()));

  compute_curvature_ = new QAction(tr("Compute Curvature"), this);
  connect(compute_curvature_, SIGNAL(triggered()),
          this, SLOT(ComputeCurvature()));

  compute_snake_length_ = new QAction(tr("Compute Snake Length"), this);
  connect(compute_snake_length_, SIGNAL(triggered()),
          this, SLOT(ComputeSnakeLength()));

  compute_all_ = new QAction(tr("Compute All"), this);
  connect(compute_all_, SIGNAL(triggered()),
          this, SLOT(ComputeAll()));

  show_analysis_options_ = new QAction(tr("Options"), this);
  connect(show_analysis_options_, SIGNAL(triggered()),
          this, SLOT(ShowAnalysisOptions()));
}

void MainWindow::CreateToolsMenuActions() {
  show_parameters_ = new QAction(tr("&Parameters"), this);
  show_parameters_->setIcon(QIcon(":/icon/Settings.png"));
  show_parameters_->setShortcut(Qt::CTRL + Qt::Key_P);
  connect(show_parameters_, SIGNAL(triggered()),
          this, SLOT(ShowParametersDialog()));

  load_viewpoint_ = new QAction(tr("Load Viewpoint"), this);
  connect(load_viewpoint_, SIGNAL(triggered()),
          this, SLOT(LoadViewpoint()));

  save_viewpoint_ = new QAction(tr("Save Viewpoint"), this);
  connect(save_viewpoint_, SIGNAL(triggered()),
          this, SLOT(SaveViewpoint()));

  save_snapshot_ = new QAction(tr("Save Snapshot"), this);
  connect(save_snapshot_, SIGNAL(triggered()), this, SLOT(SaveSnapshot()));
}

void MainWindow::CreateHelpMenuActions() {
  about_soax_ = new QAction(tr("About SOAX"), this);
  about_soax_->setStatusTip(tr("About SOAX"));
  connect(about_soax_, SIGNAL(triggered()), this, SLOT(AboutSOAX()));

  about_qt_ = new QAction(tr("About Qt"), this);
  about_qt_->setStatusTip(tr("About Qt"));
  connect(about_qt_, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void MainWindow::CreateMenus() {
  file_ = menuBar()->addMenu(tr("&File"));
  file_->addAction(open_image_);
  file_->addAction(save_as_isotropic_image_);
  file_->addAction(load_parameters_);
  file_->addAction(save_parameters_);
  file_->addAction(load_snakes_);
  file_->addAction(save_snakes_);
  file_->addAction(load_jfilament_snakes_);
  file_->addAction(save_jfilament_snakes_);
  file_->addAction(compare_snakes_);
  file_->addAction(compare_another_snakes_);
  file_->addAction(close_session_);
  file_->addAction(exit_);

  edit_ = menuBar()->addMenu(tr("&Edit"));
  snake_edit_submenu_ = edit_->addMenu(tr("Edit Mode"));
  snake_edit_submenu_->addAction(toggle_none_);
  snake_edit_submenu_->addAction(toggle_delete_snake_);
  snake_edit_submenu_->addAction(toggle_trim_tip_);
  snake_edit_submenu_->addAction(toggle_extend_tip_);
  snake_edit_submenu_->addAction(toggle_trim_body_);
  snake_edit_submenu_->addAction(toggle_delete_junction_);
  edit_->addAction(edit_snake_);

  view_ = menuBar()->addMenu(tr("&View"));
  view_->addAction(toggle_planes_);
  view_->addAction(toggle_mip_);
  view_->addAction(toggle_orientation_marker_);
  view_->addAction(toggle_corner_text_);
  view_->addAction(toggle_bounding_box_);
  view_->addAction(toggle_cube_axes_);
  view_->addAction(toggle_snakes_);
  view_->addAction(toggle_junctions_);
  view_->addAction(toggle_clip_);
  view_->addAction(toggle_color_azimuthal_);
  view_->addAction(toggle_color_polar_);
  view_->addAction(show_view_options_);

  process_ = menuBar()->addMenu(tr("&Process"));
  process_->addAction(initialize_snakes_);
  process_->addAction(deform_snakes_);
  process_->addAction(deform_snakes_in_action_);
  process_->addAction(deform_one_snake_);
  process_->addAction(cut_snakes_);
  process_->addAction(group_snakes_);

  analysis_ = menuBar()->addMenu(tr("&Analysis"));
  analysis_->addAction(compute_spherical_orientation_);
  analysis_->addAction(compute_radial_orientation_);
  analysis_->addAction(compute_point_density_);
  analysis_->addAction(compute_curvature_);
  analysis_->addAction(compute_snake_length_);
  analysis_->addAction(compute_all_);
  analysis_->addAction(show_analysis_options_);

  tools_ = menuBar()->addMenu(tr("&Tools"));
  tools_->addAction(show_parameters_);
  tools_->addAction(load_viewpoint_);
  tools_->addAction(save_viewpoint_);
  tools_->addAction(save_snapshot_);

  help_ = menuBar()->addMenu(tr("&Help"));
  help_->addAction(about_soax_);
  help_->addAction(about_qt_);
}

void MainWindow::CreateToolBar() {
  toolbar_ = addToolBar(tr("shortcuts"));
  toolbar_->addAction(open_image_);
  toolbar_->addAction(load_parameters_);
  toolbar_->addAction(save_snakes_);

  toolbar_->addSeparator();
  toolbar_->addAction(toggle_planes_);
  toolbar_->addAction(toggle_mip_);

  toolbar_->addSeparator();
  toolbar_->addAction(toggle_snakes_);
  toolbar_->addAction(toggle_junctions_);
  toolbar_->addAction(toggle_clip_);

  toolbar_->addSeparator();
  toolbar_->addAction(initialize_snakes_);
  toolbar_->addAction(deform_snakes_);
  toolbar_->addAction(cut_snakes_);
  toolbar_->addAction(group_snakes_);
  toolbar_->addAction(deform_one_snake_);

  toolbar_->addSeparator();
  toolbar_->addAction(toggle_none_);
  toolbar_->addAction(toggle_delete_snake_);
  toolbar_->addAction(toggle_trim_tip_);
  toolbar_->addAction(toggle_extend_tip_);
  toolbar_->addAction(toggle_trim_body_);
  toolbar_->addAction(toggle_delete_junction_);

  toolbar_->addSeparator();
  toolbar_->addAction(show_parameters_);

  toolbar_->addSeparator();
  toolbar_->addAction(close_session_);
}

void MainWindow::CreateStatusBar() {
  progress_bar_ = new QProgressBar;
  statusBar()->addPermanentWidget(progress_bar_);
}

void MainWindow::ResetActions() {
  save_as_isotropic_image_->setEnabled(false);
  load_snakes_->setEnabled(false);
  save_snakes_->setEnabled(false);
  load_jfilament_snakes_->setEnabled(false);
  save_jfilament_snakes_->setEnabled(false);
  compare_snakes_->setEnabled(false);
  compare_another_snakes_->setEnabled(false);
  close_session_->setEnabled(false);

  toggle_delete_snake_->setEnabled(false);
  toggle_trim_tip_->setEnabled(false);
  toggle_extend_tip_->setEnabled(false);
  toggle_trim_body_->setEnabled(false);
  toggle_delete_junction_->setEnabled(false);
  edit_snake_->setEnabled(false);

  toggle_planes_->setEnabled(false);
  toggle_mip_->setEnabled(false);
  toggle_orientation_marker_->setEnabled(false);
  toggle_corner_text_->setEnabled(false);
  toggle_bounding_box_->setEnabled(false);
  toggle_cube_axes_->setEnabled(false);
  toggle_snakes_->setEnabled(false);
  toggle_junctions_->setEnabled(false);
  toggle_clip_->setEnabled(false);
  toggle_color_azimuthal_->setEnabled(false);
  toggle_color_polar_->setEnabled(false);
  show_view_options_->setEnabled(false);

  initialize_snakes_->setEnabled(false);
  deform_snakes_->setEnabled(false);
  deform_snakes_in_action_->setEnabled(false);
  deform_one_snake_->setEnabled(false);
  cut_snakes_->setEnabled(false);
  group_snakes_->setEnabled(false);

  compute_spherical_orientation_->setEnabled(false);
  compute_radial_orientation_->setEnabled(false);
  compute_point_density_->setEnabled(false);
  compute_curvature_->setEnabled(false);
  compute_snake_length_->setEnabled(false);
  compute_all_->setEnabled(false);
  show_analysis_options_->setEnabled(false);

  show_parameters_->setEnabled(false);
  load_viewpoint_->setEnabled(false);
  save_viewpoint_->setEnabled(false);
  save_snapshot_->setEnabled(false);
}

void MainWindow::OpenImage() {
  QString dir = this->GetLastDirectory(image_filename_);
  QString filename = QFileDialog::getOpenFileName(
      this, tr("Open Image"), dir,
      tr("Image Files (*.tif *.tiff *.mhd *.mha *.png *.jpg *.bmp)"));

  if (filename.isEmpty()) return;
  image_filename_ = filename.toStdString();
  this->setWindowTitle(filename.prepend("SOAX - "));

  multisnake_->LoadImage(image_filename_);
  viewer_->SetupImage(multisnake_->image());
  toggle_planes_->setChecked(true);
  if (multisnake_->dim() > 2)
    toggle_mip_->setChecked(true);
  toggle_orientation_marker_->setChecked(true);
  toggle_corner_text_->setChecked(true);
  toggle_bounding_box_->setChecked(false);
  toggle_cube_axes_->setChecked(false);
  statusBar()->showMessage(tr("Image loaded."), message_timeout_);

  view_options_dialog_->SetWindow(viewer_->window());
  view_options_dialog_->SetLevel(viewer_->level());
  view_options_dialog_->SetMinIntensity(viewer_->mip_min_intensity());
  view_options_dialog_->SetMaxIntensity(viewer_->mip_max_intensity());
  view_options_dialog_->SetClipSpan(viewer_->clip_span());
  view_options_dialog_->SetColorSegmentStep(viewer_->color_segment_step());

  PointType center = multisnake_->GetImageCenter();
  analysis_options_dialog_->SetImageCenter(center);
  DataContainer center_coordinates;
  for (int i = 0; i < kDimension; i++) {
    center_coordinates.push_back(center[i]);
  }
  // double radius = Minimum(center_coordinates);
  double radius = multisnake_->GetImageDiagonal() / 2.0;
  // std::cout << radius <<  std::endl;
  analysis_options_dialog_->SetRadius(radius);

  open_image_->setEnabled(false);
  save_as_isotropic_image_->setEnabled(true);
  load_snakes_->setEnabled(true);
  load_jfilament_snakes_->setEnabled(true);
  compare_snakes_->setEnabled(true);
  close_session_->setEnabled(true);
  toggle_planes_->setEnabled(true);
  toggle_mip_->setEnabled(true);
  toggle_orientation_marker_->setEnabled(true);
  toggle_corner_text_->setEnabled(true);
  toggle_bounding_box_->setEnabled(true);
  toggle_cube_axes_->setEnabled(true);
  show_view_options_->setEnabled(true);
  initialize_snakes_->setEnabled(true);
  show_parameters_->setEnabled(true);
  load_viewpoint_->setEnabled(true);
  save_viewpoint_->setEnabled(true);
  save_snapshot_->setEnabled(true);
}

void MainWindow::SaveAsIsotropicImage() {
  bool ok;
  double z_spacing = QInputDialog::getDouble(
      this, tr("Set Z spacing"), tr("Z Spacing (relative to X/Y)"),
      1.0, 0.1, 10, 4, &ok);
  if (!ok) return;

  if (std::fabs(z_spacing - 1) > kEpsilon) {
    QString dir = this->GetLastDirectory(image_filename_);
    QString filename = QFileDialog::getSaveFileName(
        this, tr("Save as Isotropic Image"), dir);
    if (filename.isEmpty()) return;
    image_filename_ = filename.toStdString();

    multisnake_->SaveAsIsotropicImage(image_filename_, z_spacing);
    statusBar()->showMessage(tr("Image has been resampled and saved."),
                             message_timeout_);
  } else {
    statusBar()->showMessage(tr("Image is already isotropic! Abort."),
                             message_timeout_);
  }
}

void MainWindow::LoadParameters() {
  QString dir = this->GetLastDirectory(parameter_filename_);
  QString filename = QFileDialog::getOpenFileName(
      this, tr("Open Parameter File"), dir, tr("Text Files (*.txt)"));
  if (filename.isEmpty()) return;
  parameter_filename_ = filename.toStdString();

  multisnake_->LoadParameters(parameter_filename_);
  statusBar()->showMessage(tr("Parameters loaded."), message_timeout_);
}

void MainWindow::SaveParameters() {
  QString dir = this->GetLastDirectory(parameter_filename_);
  QString filename = QFileDialog::getSaveFileName(
      this, tr("Save Parameters"), dir, tr("Text Files (*.txt)"));
  if (filename.isEmpty()) return;
  parameter_filename_ = filename.toStdString();

  multisnake_->SaveParameters(parameter_filename_);
  statusBar()->showMessage("Parameters saved.", message_timeout_);
}

void MainWindow::LoadSnakes() {
  QString dir = this->GetLastDirectory(snake_filename_);
  QString filename = QFileDialog::getOpenFileName(
      this, tr("Load Snakes"), dir, tr("Text Files (*.txt)"));
  if (filename.isEmpty()) return;
  snake_filename_ = filename.toStdString();

  multisnake_->LoadConvergedSnakes(snake_filename_);
  unsigned num_snakes = multisnake_->GetNumberOfConvergedSnakes();
  QString msg = "Number of snakes loaded: " + QString::number(num_snakes);
  statusBar()->showMessage(msg, message_timeout_);

  viewer_->RemoveJunctions();
  viewer_->RemoveSnakes();
  viewer_->SetupSnakes(multisnake_->converged_snakes());
  viewer_->SetupSnakes(multisnake_->comparing_snakes1(), 1);
  viewer_->SetupSnakes(multisnake_->comparing_snakes2(), 2);
  toggle_snakes_->setChecked(true);
  viewer_->set_snake_filename(snake_filename_);
  viewer_->SetupUpperRightCornerText();
  toggle_corner_text_->setChecked(true);
  viewer_->SetupJunctions(multisnake_->GetJunctions());
  toggle_junctions_->setChecked(true);
  viewer_->Render();

  save_snakes_->setEnabled(true);
  save_jfilament_snakes_->setEnabled(true);
  compare_snakes_->setEnabled(true);
  toggle_delete_snake_->setEnabled(true);
  toggle_trim_tip_->setEnabled(true);
  toggle_extend_tip_->setEnabled(true);
  toggle_trim_body_->setEnabled(true);
  toggle_delete_junction_->setEnabled(true);
  edit_snake_->setEnabled(true);
  initialize_snakes_->setEnabled(false);
  deform_one_snake_->setEnabled(false);
  toggle_snakes_->setEnabled(true);
  toggle_junctions_->setEnabled(true);
  toggle_clip_->setEnabled(true);
  toggle_color_azimuthal_->setEnabled(true);
  toggle_color_polar_->setEnabled(true);
  compute_spherical_orientation_->setEnabled(true);
  compute_radial_orientation_->setEnabled(true);
  compute_point_density_->setEnabled(true);
  compute_curvature_->setEnabled(true);
  compute_snake_length_->setEnabled(true);
  compute_all_->setEnabled(true);
  show_analysis_options_->setEnabled(true);
}

void MainWindow::SaveSnakes() {
  QString dir = this->GetLastDirectory(snake_filename_);
  QString filename = QFileDialog::getSaveFileName(
      this, tr("Save Snakes"), dir, tr("Text Files (*.txt)"));
  if (filename.isEmpty()) return;
  snake_filename_ = filename.toStdString();

  multisnake_->SaveSnakes(multisnake_->converged_snakes(),
                          snake_filename_);
  statusBar()->showMessage(tr("Snakes are saved."), message_timeout_);
}

void MainWindow::LoadJFilamentSnakes() {
  QString dir = this->GetLastDirectory(snake_filename_);
  QString filename = QFileDialog::getOpenFileName(
      this, tr("Load JFilament Snakes"), dir, tr("Text Files (*.txt)"));
  if (filename.isEmpty()) return;
  snake_filename_ = filename.toStdString();

  multisnake_->LoadGroundTruthSnakes(snake_filename_);
  unsigned num_snakes = multisnake_->GetNumberOfComparingSnakes1();
  QString msg = "# JFilament snakes loaded: " +
      QString::number(num_snakes);
  statusBar()->showMessage(msg, message_timeout_);
  viewer_->RemoveSnakes();
  viewer_->SetupSnakes(multisnake_->converged_snakes());

  viewer_->SetupSnakes(multisnake_->comparing_snakes1(), 1);
  viewer_->SetupSnakes(multisnake_->comparing_snakes2(), 2);
  toggle_snakes_->setChecked(true);
  viewer_->set_comapring_snake_filename1(snake_filename_);
  viewer_->SetupUpperRightCornerText();
  toggle_corner_text_->setChecked(true);
  viewer_->Render();

  save_snakes_->setEnabled(true);
  save_jfilament_snakes_->setEnabled(true);
  compare_snakes_->setEnabled(true);
  toggle_delete_snake_->setEnabled(true);
  toggle_trim_tip_->setEnabled(true);
  toggle_extend_tip_->setEnabled(true);
  toggle_trim_body_->setEnabled(true);
  toggle_delete_junction_->setEnabled(true);
  edit_snake_->setEnabled(true);
  toggle_snakes_->setEnabled(true);
  toggle_clip_->setEnabled(true);
  toggle_color_azimuthal_->setEnabled(true);
  toggle_color_polar_->setEnabled(true);
  compute_spherical_orientation_->setEnabled(true);
  compute_radial_orientation_->setEnabled(true);
  compute_point_density_->setEnabled(true);
  compute_curvature_->setEnabled(true);
  compute_all_->setEnabled(true);
  show_analysis_options_->setEnabled(true);
}

void MainWindow::SaveJFilamentSnakes() {
  QString dir = this->GetLastDirectory(snake_filename_);
  QString filename = QFileDialog::getSaveFileName(
      this, tr("Save JFilament Snakes"), dir, tr("Text Files (*.txt)"));
  if (filename.isEmpty()) return;
  snake_filename_ = filename.toStdString();

  multisnake_->SaveConvergedSnakesAsJFilamentFormat(snake_filename_);
  statusBar()->showMessage(tr("Snakes are saved in JFilament format"),
                           message_timeout_);
}

void MainWindow::CompareSnakes() {
  QString dir = this->GetLastDirectory(snake_filename_);
  QString filename = QFileDialog::getOpenFileName(
      this, tr("Compare Snakes"), dir, tr("Text Files (*.txt)"));
  if (filename.isEmpty()) return;
  snake_filename_ = filename.toStdString();

  multisnake_->LoadComparingSnakes1(snake_filename_);
  unsigned num_snakes = multisnake_->GetNumberOfComparingSnakes1();
  QString msg = "Number of comparing snakes loaded: " +
      QString::number(num_snakes);
  statusBar()->showMessage(msg, message_timeout_);
  viewer_->RemoveSnakes();
  viewer_->SetupSnakes(multisnake_->converged_snakes());
  viewer_->SetupSnakes(multisnake_->comparing_snakes1(), 1);
  viewer_->SetupSnakes(multisnake_->comparing_snakes2(), 2);
  toggle_snakes_->setChecked(true);
  viewer_->set_comapring_snake_filename1(snake_filename_);
  viewer_->SetupUpperRightCornerText();
  toggle_corner_text_->setChecked(true);
  viewer_->Render();

  compare_another_snakes_->setEnabled(true);
  toggle_delete_snake_->setEnabled(false);
  toggle_trim_tip_->setEnabled(false);
  toggle_extend_tip_->setEnabled(false);
  toggle_trim_body_->setEnabled(false);
  toggle_delete_junction_->setEnabled(false);
  edit_snake_->setEnabled(false);
  toggle_snakes_->setEnabled(true);
  toggle_clip_->setEnabled(true);
  toggle_color_azimuthal_->setEnabled(false);
  toggle_color_polar_->setEnabled(false);
}

void MainWindow::CompareAnotherSnakes() {
  QString dir = this->GetLastDirectory(snake_filename_);
  QString filename = QFileDialog::getOpenFileName(
      this, tr("Compare Other Snakes"), dir, tr("Text Files (*.txt)"));
  if (filename.isEmpty()) return;
  snake_filename_ = filename.toStdString();

  multisnake_->LoadComparingSnakes2(snake_filename_);
  unsigned num_snakes = multisnake_->GetNumberOfComparingSnakes2();
  QString msg = "Number of other comparing snakes loaded: " +
      QString::number(num_snakes);
  statusBar()->showMessage(msg, message_timeout_);
  viewer_->RemoveSnakes();
  viewer_->SetupSnakes(multisnake_->converged_snakes());
  viewer_->SetupSnakes(multisnake_->comparing_snakes1(), 1);
  viewer_->SetupSnakes(multisnake_->comparing_snakes2(), 2);
  toggle_snakes_->setChecked(true);
  viewer_->set_comapring_snake_filename2(snake_filename_);
  viewer_->SetupUpperRightCornerText();
  toggle_corner_text_->setChecked(true);
  viewer_->Render();

  toggle_snakes_->setEnabled(true);
}

void MainWindow::CloseSession() {
  parameters_dialog_->SetCurrentParameters(multisnake_);

  toggle_planes_->setChecked(false);
  toggle_mip_->setChecked(false);
  toggle_orientation_marker_->setChecked(false);
  toggle_corner_text_->setChecked(false);
  toggle_bounding_box_->setChecked(false);;
  toggle_cube_axes_->setChecked(false);
  toggle_snakes_->setChecked(false);
  toggle_junctions_->setChecked(false);
  toggle_clip_->setChecked(false);
  toggle_color_azimuthal_->setChecked(false);
  toggle_color_polar_->setChecked(false);

  viewer_->Reset();
  viewer_->Render();
  multisnake_->Reset();

  this->ResetActions();
  open_image_->setEnabled(true);
  setWindowTitle("SOAX");
}

void MainWindow::EditSnake() {
  if (toggle_delete_snake_->isChecked()) {
    multisnake_->DeleteSnakes(viewer_->selected_snakes());
    viewer_->RemoveSelectedSnakes();
  } else if (toggle_trim_tip_->isChecked()) {
    viewer_->TrimTip();
  } else if (toggle_extend_tip_->isChecked()) {
    viewer_->ExtendTip();
  } else if (toggle_trim_body_->isChecked()) {
    viewer_->TrimBody();
  } else if (toggle_delete_junction_->isChecked()) {
    viewer_->RemoveSelectedJunctions(multisnake_->junctions());
  }
  viewer_->Render();
  deform_one_snake_->setEnabled(true);
}

void MainWindow::ShowViewOptions() {
  if (view_options_dialog_->exec()) {
    double window = view_options_dialog_->GetWindow();
    double level = view_options_dialog_->GetLevel();
    viewer_->UpdateWindowLevel(window, level);
    double min_intensity = view_options_dialog_->GetMinIntensity();
    double max_intensity = view_options_dialog_->GetMaxIntensity();
    viewer_->UpdateMIPIntensityRange(min_intensity, max_intensity);
    viewer_->set_clip_span(view_options_dialog_->GetClipSpan());

    if (viewer_->color_segment_step() !=
        view_options_dialog_->GetColorSegmentStep()) {
      viewer_->set_color_segment_step(
          view_options_dialog_->GetColorSegmentStep());
      if (toggle_color_azimuthal_->isChecked()) {
        viewer_->ColorByAzimuthalAngle(true);
      } else if (toggle_color_polar_->isChecked()) {
        viewer_->ColorByPolarAngle(true);
      }
    }
    view_options_dialog_->DisableOKButton();
  }
  viewer_->Render();
}

void MainWindow::InitializeSnakes() {
  multisnake_->ComputeImageGradient();

  multisnake_->InitializeSnakes();
  QString msg = QString::number(multisnake_->GetNumberOfInitialSnakes()) +
      " snakes initialized.";
  statusBar()->showMessage(msg, message_timeout_);

  viewer_->RemoveSnakes();
  viewer_->SetupSnakes(multisnake_->initial_snakes());
  toggle_snakes_->setChecked(true);
  viewer_->Render();

  save_snakes_->setEnabled(true);
  save_jfilament_snakes_->setEnabled(true);
  toggle_snakes_->setEnabled(true);
  toggle_clip_->setEnabled(true);
  deform_snakes_->setEnabled(true);
  deform_snakes_in_action_->setEnabled(true);
}

void MainWindow::DeformSnakes() {
  std::cout << "============ Parameters ============" << std::endl;
  multisnake_->WriteParameters(std::cout);
  std::cout << "====================================" << std::endl;

  progress_bar_->setMaximum(multisnake_->GetNumberOfInitialSnakes());
  connect(multisnake_, SIGNAL(ExtractionProgressed(int)),
          progress_bar_, SLOT(setValue(int)));

  time_t start, end;
  time(&start);
  multisnake_->DeformSnakes();
  time(&end);
  double time_elasped = difftime(end, start) / 60.0;
  std::cout << "Extraction complete. Evolution time: "
            << time_elasped << "m" << std::endl;
  unsigned num_snakes = multisnake_->GetNumberOfConvergedSnakes();
  QString msg = "# converged snakes: " + QString::number(num_snakes);
  statusBar()->showMessage(msg, message_timeout_);

  viewer_->RemoveSnakes();
  viewer_->SetupSnakes(multisnake_->converged_snakes());
  toggle_snakes_->setChecked(true);
  viewer_->Render();

  initialize_snakes_->setEnabled(false);
  deform_snakes_->setEnabled(false);
  save_snakes_->setEnabled(true);
  save_jfilament_snakes_->setEnabled(true);
  compare_snakes_->setEnabled(true);
  toggle_snakes_->setEnabled(true);
  toggle_clip_->setEnabled(true);
  toggle_color_azimuthal_->setEnabled(true);
  toggle_color_polar_->setEnabled(true);

  toggle_delete_snake_->setEnabled(true);
  toggle_trim_tip_->setEnabled(true);
  toggle_extend_tip_->setEnabled(true);
  toggle_trim_body_->setEnabled(true);
  edit_snake_->setEnabled(true);

  cut_snakes_->setEnabled(true);
  compute_spherical_orientation_->setEnabled(true);
  compute_radial_orientation_->setEnabled(true);
  compute_point_density_->setEnabled(true);
  compute_curvature_->setEnabled(true);
  compute_snake_length_->setEnabled(true);
  compute_all_->setEnabled(true);
  show_analysis_options_->setEnabled(true);
}

void MainWindow::DeformSnakesInAction() {
  viewer_->RemoveSnakes();
  progress_bar_->setMaximum(multisnake_->GetNumberOfInitialSnakes());
  unsigned ncompleted = 0;
  
  // clear converged_snakes_grid_ while keeping the size the same
  multisnake_->ClearConvergedSnakesGrid();
  
  while (!multisnake_->initial_snakes().empty()) {
    Snake *s = multisnake_->PopLastInitialSnake();
    viewer_->SetupSnake(s, 0);
    viewer_->Render();
    s->Evolve(multisnake_->solver_bank(), multisnake_->converged_snakes(),
              Snake::iterations_per_press(), multisnake_->dim(), multisnake_->converged_snakes_grid());

    if (s->viable()) {
      if (s->converged()) {
        viewer_->SetupSnake(s, 0);
        viewer_->ChangeSnakeColor(s, Viewer::Yellow());
        multisnake_->AddConvergedSnake(s);
        
        for(int iterate_snakes = 0; iterate_snakes < s->GetSize(); iterate_snakes++) {
           double curr_x_val = s->GetX(iterate_snakes);
           double curr_y_val = s->GetY(iterate_snakes);
           
           int org_x_grid = (int)(curr_x_val / Snake::overlap_threshold());
           int org_y_grid = (int)(curr_y_val / Snake::overlap_threshold());
          
           multisnake_->AddConvergedSnakeIndexesToGrid(org_x_grid, org_y_grid, multisnake_->GetNumberOfConvergedSnakes()-1, iterate_snakes);
        }
        ncompleted++;
      } else {
        multisnake_->AddInitialSnake(s);
      }
    } else {
      multisnake_->AddSubsnakesToInitialSnakes(s);
      viewer_->RemoveSnake(s);
      delete s;
      ncompleted++;
    }
    progress_bar_->setValue(ncompleted);
    qApp->processEvents();
    viewer_->Render();
    std::cout << "\rRemaining: " << std::setw(6)
              << multisnake_->GetNumberOfInitialSnakes() << std::flush;
  }
  statusBar()->showMessage(tr("Evolution complete."));
  viewer_->RemoveSnakes();
  viewer_->SetupSnakes(multisnake_->converged_snakes());
  toggle_snakes_->setChecked(true);
  viewer_->Render();

  initialize_snakes_->setEnabled(false);
  deform_snakes_->setEnabled(false);
  deform_snakes_in_action_->setEnabled(false);

  save_snakes_->setEnabled(true);
  save_jfilament_snakes_->setEnabled(true);
  compare_snakes_->setEnabled(true);
  toggle_delete_snake_->setEnabled(true);
  toggle_trim_tip_->setEnabled(true);
  toggle_extend_tip_->setEnabled(true);
  toggle_trim_body_->setEnabled(true);
  edit_snake_->setEnabled(true);

  toggle_snakes_->setEnabled(true);
  toggle_clip_->setEnabled(true);
  toggle_color_azimuthal_->setEnabled(true);
  toggle_color_polar_->setEnabled(true);
  cut_snakes_->setEnabled(true);
  compute_spherical_orientation_->setEnabled(true);
  compute_radial_orientation_->setEnabled(true);
  compute_point_density_->setEnabled(true);
  compute_curvature_->setEnabled(true);
  compute_all_->setEnabled(true);
  show_analysis_options_->setEnabled(true);
}

void MainWindow::DeformOneSnake() {
  if (viewer_->trimmed_snake()) {
    if (!multisnake_->external_force()) {
      multisnake_->ComputeImageGradient();
    }

    viewer_->trimmed_snake()->EvolveWithTipFixed(
        multisnake_->solver_bank(), Snake::iterations_per_press(),
        multisnake_->dim());

    if (viewer_->trimmed_snake()->converged()) {
      statusBar()->showMessage(tr("Snake is converged."));
    }
    viewer_->SetupSnake(viewer_->trimmed_snake(), 0);
    viewer_->Render();
  } else {
    std::cerr << "No edited snake!" << std::endl;
  }
}

void MainWindow::CutSnakes() {
  multisnake_->CutSnakesAtTJunctions();
  viewer_->RemoveSnakes();
  viewer_->SetupSnakes(multisnake_->converged_snakes());
  toggle_snakes_->setChecked(true);
  viewer_->Render();
  statusBar()->showMessage(tr("Snakes are cut at junctions."));
  cut_snakes_->setEnabled(false);
  group_snakes_->setEnabled(true);
  save_snakes_->setEnabled(true);
  save_jfilament_snakes_->setEnabled(true);

  compute_spherical_orientation_->setEnabled(true);
  compute_radial_orientation_->setEnabled(true);
  compute_point_density_->setEnabled(true);
  compute_curvature_->setEnabled(true);
  compute_all_->setEnabled(true);
  show_analysis_options_->setEnabled(true);
}

void MainWindow::GroupSnakes() {
  multisnake_->GroupSnakes();
  viewer_->RemoveSnakes();
  viewer_->SetupSnakes(multisnake_->converged_snakes());
  toggle_snakes_->setChecked(true);
  viewer_->SetupJunctions(multisnake_->GetJunctions());
  toggle_junctions_->setChecked(true);
  viewer_->Render();
  statusBar()->showMessage(tr("Snakes are reconfigured."));

  group_snakes_->setEnabled(false);
  toggle_junctions_->setEnabled(true);
  save_snakes_->setEnabled(true);
  save_jfilament_snakes_->setEnabled(true);
  compare_snakes_->setEnabled(true);
  toggle_delete_junction_->setEnabled(true);

  compute_spherical_orientation_->setEnabled(true);
  compute_radial_orientation_->setEnabled(true);
  compute_point_density_->setEnabled(true);
  compute_curvature_->setEnabled(true);
  compute_all_->setEnabled(true);
  show_analysis_options_->setEnabled(true);
}

void MainWindow::ComputeSphericalOrientation() {
  QString dir = this->GetLastDirectory(analysis_filename_);
  QString filename = QFileDialog::getSaveFileName(
      this, tr("Save spherical orientation"), dir, tr("CSV files (*.csv)"));
  if (filename.isEmpty()) return;
  analysis_filename_ = filename.toStdString();
  if (WriteSphericalOrientation(analysis_filename_))
    statusBar()->showMessage(tr("Spherical orientation file saved."));
}

bool MainWindow::WriteSphericalOrientation(const std::string &filename) {
  PointType center;
  if (!analysis_options_dialog_->GetImageCenter(&center)) {
    ShowErrorDialog("Image center is invalid!");
    return false;
  }

  double inside_percentage = 1.0;
  if (!analysis_options_dialog_->GetInsideRatio(&inside_percentage)) {
    ShowErrorDialog("Inside ratio is invalid!");
    return false;
  }

  double radius = 0.0;
  if (!analysis_options_dialog_->GetRadius(&radius)) {
    ShowErrorDialog("Radius is invalid!");
    return false;
  }

  double padding = 0.0;
  if (analysis_options_dialog_->ExcludeBoundaryChecked()) {
    padding = 2.0;
  }

  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    ShowErrorDialog("Open file failed!");
    return false;
  }

  double max_r = inside_percentage * radius;
  multisnake_->ComputeSphericalOrientation(center, max_r, padding, outfile);
  outfile.close();
  return true;
}

void MainWindow::ComputeRadialOrientation() {
  QString dir = this->GetLastDirectory(analysis_filename_);
  QString filename = QFileDialog::getSaveFileName(
      this, tr("Save radial orientation"), dir, tr("CSV files (*.csv)"));
  if (filename.isEmpty()) return;
  analysis_filename_ = filename.toStdString();
  if (WriteRadialOrientation(analysis_filename_))
    statusBar()->showMessage(tr("Radial orientation file saved."));
}

bool MainWindow::WriteRadialOrientation(const std::string &filename) {
  PointType center;
  if (!analysis_options_dialog_->GetImageCenter(&center)) {
    ShowErrorDialog("Image center is invalid!");
    return false;
  }
  double pixel_size = 1.0;
  if (!analysis_options_dialog_->GetPixelSize(&pixel_size)) {
    ShowErrorDialog("Pixel size is invalid!");
    return false;
  }

  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    ShowErrorDialog("Open file failed!");
    return false;
  }

  multisnake_->ComputeRadialOrientation(center, pixel_size, outfile);
  outfile.close();
  return true;
}

void MainWindow::ComputePointDensity() {
  QString dir = this->GetLastDirectory(analysis_filename_);
  QString filename = QFileDialog::getSaveFileName(
      this, tr("Save SOAC point density/intensity"), dir,
      tr("CSV files (*.csv)"));
  if (filename.isEmpty()) return;
  analysis_filename_ = filename.toStdString();
  if (WritePointDensity(analysis_filename_))
    statusBar()->showMessage(tr("SOAC point density/intensity file saved."));
}

bool MainWindow::WritePointDensity(const std::string &filename) {
  PointType center;
  if (!analysis_options_dialog_->GetImageCenter(&center)) {
    ShowErrorDialog("Image center is invalid!");
    return false;
  }

  double radius = 0.0;
  if (!analysis_options_dialog_->GetRadius(&radius)) {
    ShowErrorDialog("Radius is invalid!");
    return false;
  }

  double inside_percentage = 1.0;
  if (!analysis_options_dialog_->GetInsideRatio(&inside_percentage)) {
    ShowErrorDialog("Inside ratio is invalid!");
    return false;
  }

  // double inside_percentage = 1.1;
  double max_r = inside_percentage * radius;
  double pixel_size = 1.0;
  if (!analysis_options_dialog_->GetPixelSize(&pixel_size)) {
    ShowErrorDialog("Pixel size is invalid!");
    return false;
  }

  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    ShowErrorDialog("Open file failed!");
    return false;
  }

  multisnake_->ComputePointDensityAndIntensity(center, max_r,
                                               pixel_size, outfile);
  outfile.close();
  return true;
}

void MainWindow::ComputeCurvature() {
  QString dir = this->GetLastDirectory(analysis_filename_);
  QString filename = QFileDialog::getSaveFileName(
      this, tr("Save curvature"), dir, tr("CSV files (*.csv)"));
  if (filename.isEmpty()) return;
  analysis_filename_ = filename.toStdString();
  if (WriteCurvature(analysis_filename_))
    statusBar()->showMessage(tr("Curvature file saved."));
}

bool MainWindow::WriteCurvature(const std::string &filename) {
  int coarse_graining = 8;
  if (!analysis_options_dialog_->GetCoarseGraining(&coarse_graining)) {
    ShowErrorDialog("Coarse graining is invalid!");
    return false;
  }

  double pixel_size = 1.0;
  if (!analysis_options_dialog_->GetPixelSize(&pixel_size)) {
    ShowErrorDialog("Pixel size is invalid!");
    return false;
  }

  double padding = 0.0;
  if (analysis_options_dialog_->ExcludeBoundaryChecked()) {
    padding = 2.0;
  }

  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    ShowErrorDialog("Open file failed!");
    return false;
  }

  multisnake_->ComputeCurvature(coarse_graining, pixel_size,
                                padding, outfile);
  outfile.close();
  return true;
}

void MainWindow::ComputeSnakeLength() {
  QString dir = this->GetLastDirectory(analysis_filename_);
  QString filename = QFileDialog::getSaveFileName(
      this, tr("Save Snake Length"), dir, tr("CSV files (*.csv)"));
  if (filename.isEmpty()) return;
  analysis_filename_ = filename.toStdString();
  if (WriteSnakeLength(analysis_filename_))
    statusBar()->showMessage(tr("Snake length file saved."));
}

bool MainWindow::WriteSnakeLength(const std::string &filename) {
  double pixel_size = 1.0;
  if (!analysis_options_dialog_->GetPixelSize(&pixel_size)) {
    ShowErrorDialog("Pixel size is invalid!");
    return false;
  }

  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    ShowErrorDialog("Open file failed!");
    return false;
  }
  outfile << "Index,Length(um)" << std::endl;
  multisnake_->ComputeSnakeLength(pixel_size, outfile);
  outfile.close();
  return true;
}

void MainWindow::ComputeAll() {
  QString prev_dir = this->GetLastDirectory(analysis_filename_);
  QString curr_dir = QFileDialog::getExistingDirectory(
      this, tr("Open an Analysis Folder"), prev_dir,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
  analysis_filename_ = curr_dir.toStdString();

  std::string spherical_orientation_file_path =
      GenerateAnalysisFilePath(curr_dir, "_spherical_orientation.csv");
  if (!WriteSphericalOrientation(spherical_orientation_file_path)) {
    ShowErrorDialog("Write spherical orientation failed!");
    return;
  }

  std::string radial_orientation_file_path =
      GenerateAnalysisFilePath(curr_dir, "_radial_orientation.csv");
  if (!WriteRadialOrientation(radial_orientation_file_path)) {
    ShowErrorDialog("Write radial orientation failed!");
    return;
  }

  std::string point_density_file_path =
      GenerateAnalysisFilePath(curr_dir, "_point_density.csv");
  if (!WritePointDensity(point_density_file_path)) {
    ShowErrorDialog("Write point density failed!");
    return;
  }

  std::string curvature_file_path =
      GenerateAnalysisFilePath(curr_dir, "_curvature.csv");
  if (!WriteCurvature(curvature_file_path)) {
    ShowErrorDialog("Write curvature failed!");
    return;
  }

  std::string length_file_path =
      GenerateAnalysisFilePath(curr_dir, "_length.csv");
  if (!WriteSnakeLength(length_file_path)) {
    ShowErrorDialog("Write snake length failed!");
    return;
  }

  // std::cout << spherical_orientation_file_path << std::endl;
  // std::cout << radial_orientation_file_path << std::endl;
  // std::cout << point_density_file_path << std::endl;
  // std::cout << curvature_file_path << std::endl;
  // std::cout << length_file_path << std::endl;

  QString msg = QString("All files are saved in ") + curr_dir;
  statusBar()->showMessage(msg);
}

std::string MainWindow::GenerateAnalysisFilePath(const QString &dir,
                                                 const std::string &str) const {
  std::string::size_type pos = image_filename_.find_last_of("/\\");
  std::string name = image_filename_.substr(pos + 1,
                                            image_filename_.size() - pos - 5) + str;
  return QDir(dir).absoluteFilePath(name.c_str()).toStdString();
}

void MainWindow::ShowAnalysisOptions() {
  if (analysis_options_dialog_->exec()) {
    analysis_options_dialog_->DisableOKButton();
  }
}

void MainWindow::ShowParametersDialog() {
  parameters_dialog_->SetCurrentParameters(multisnake_);
  if (parameters_dialog_->exec()) {
    this->SetParameters();
    parameters_dialog_->DisableOKButton();
  }
}

void MainWindow::SetParameters() {
  multisnake_->set_intensity_scaling(
      parameters_dialog_->GetIntensityScaling());
  Snake::set_intensity_scaling(multisnake_->intensity_scaling());
  multisnake_->set_sigma(parameters_dialog_->GetSigma());
  multisnake_->set_ridge_threshold(parameters_dialog_->GetRidgeThreshold());
  multisnake_->set_foreground(parameters_dialog_->GetForeground());
  multisnake_->set_background(parameters_dialog_->GetBackground());
  Snake::set_background(multisnake_->background());
  Snake::set_desired_spacing(parameters_dialog_->GetSpacing());
  multisnake_->set_initialize_z(parameters_dialog_->InitializeZ());
  Snake::set_minimum_length(parameters_dialog_->GetMinSnakeLength());
  Snake::set_max_iterations(parameters_dialog_->GetMaxIterations());
  Snake::set_change_threshold(parameters_dialog_->GetChangeThreshold());
  Snake::set_check_period(parameters_dialog_->GetCheckPeriod());
  Snake::set_iterations_per_press(
      parameters_dialog_->GetIterationsPerPress());
  multisnake_->solver_bank()->set_alpha(parameters_dialog_->GetAlpha());
  multisnake_->solver_bank()->set_beta(parameters_dialog_->GetBeta());
  multisnake_->solver_bank()->set_gamma(parameters_dialog_->GetGamma());
  Snake::set_external_factor(parameters_dialog_->GetExternalFactor());
  Snake::set_stretch_factor(parameters_dialog_->GetStretchFactor());
  Snake::set_number_of_sectors(parameters_dialog_->GetNumberOfSectors());
  Snake::set_radial_near(parameters_dialog_->GetRadialNear());
  Snake::set_radial_far(parameters_dialog_->GetRadialFar());
  Snake::set_radial_save_foreground(parameters_dialog_->GetRadialSaveForeground());
  Snake::set_z_spacing(parameters_dialog_->GetZSpacing());
  Snake::set_delta(parameters_dialog_->GetDelta());
  Snake::set_overlap_threshold(parameters_dialog_->GetOverlapThreshold());
  Snake::set_grouping_distance_threshold(
      parameters_dialog_->GetGroupingDistanceThreshold());
  Snake::set_grouping_delta(parameters_dialog_->GetGroupingDelta());
  Snake::set_direction_threshold(parameters_dialog_->GetDirectionThreshold());
  Snake::set_damp_z(parameters_dialog_->DampZ());
}

void MainWindow::LoadViewpoint() {
  QString dir = this->GetLastDirectory(viewpoint_filename_);
  QString filename = QFileDialog::getOpenFileName(
      this, tr("Load Viewpoint"), dir, tr("Text Files (*.cam *.txt)"));
  if (filename.isEmpty()) return;
  viewpoint_filename_ = filename.toStdString();

  viewer_->LoadViewpoint(viewpoint_filename_);
  statusBar()->showMessage(tr("Viewpoint loaded."), message_timeout_);
}

void MainWindow::SaveViewpoint() {
  QString dir = this->GetLastDirectory(viewpoint_filename_);
  QString filename = QFileDialog::getSaveFileName(
      this, tr("Save Viewpoint"), dir, tr("Text Files (*.cam *.txt)"));
  if (filename.isEmpty()) return;
  viewpoint_filename_ = filename.toStdString();

  viewer_->SaveViewpoint(viewpoint_filename_);
  statusBar()->showMessage("Viewpoint saved.", message_timeout_);
}

void MainWindow::SaveSnapshot() {
  QString dir = this->GetLastDirectory(snapshot_filename_);
  QString filename = QFileDialog::getSaveFileName(
      this, tr("Save Snapshot"), dir, tr("PNG Image Files (*.png)"));
  if (filename.isEmpty()) return;
  snapshot_filename_ = filename.toStdString();

  viewer_->PrintScreenAsPNGImage(snapshot_filename_);
  statusBar()->showMessage("Snapshots saved.", message_timeout_);
}

void MainWindow::AboutSOAX() {
  QMessageBox::about(
      this, tr("About SOAX"),
      tr("<h3>SOAX 3.8.0</h3>"
         "<p>Copyright &copy; Lehigh University"
         "<p>SOAX extracts curvilinear networks from 2D/3D images."
         "This work was supported by NIH grants R01GM098430 and R35GM136372."));
}

QString MainWindow::GetLastDirectory(const std::string &filename) const {
  QString dir = "..";
  if (!filename.empty()) {  // extract the last directory
    std::string::size_type pos = filename.find_last_of("/\\");
    dir = QString(filename.substr(0, pos).c_str());
  }
  return dir;
}

void MainWindow::ShowErrorDialog(const char *msg) {
  QMessageBox msg_box;
  msg_box.setText(msg);
  msg_box.setIcon(QMessageBox::Warning);
  msg_box.exec();
}

}  // namespace soax
