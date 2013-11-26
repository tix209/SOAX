#include <QtGui>
#include "QVTKWidget.h"
#include "main_window.h"
#include "multisnake.h"
#include "viewer.h"
#include "parameters_dialog.h"
#include "solver_bank.h"
#include "analysis_options_dialog.h"

namespace soax {

MainWindow::MainWindow() {
  central_widget_ = new QWidget;
  multisnake_ = new Multisnake;
  viewer_ = new Viewer;
  parameters_dialog_ = new ParametersDialog(this);
  analysis_options_dialog_ = new AnalysisOptionsDialog(this);
  QHBoxLayout* layout = new QHBoxLayout;
  layout->addWidget(viewer_->qvtk());
  central_widget_->setLayout(layout);
  setCentralWidget(central_widget_);
  setWindowIcon(QIcon(":/icon/letter-x.png"));
  setWindowTitle("SOAX");
  this->CreateActions();
  this->CreateMenus();
  this->CreateToolBar();
  this->CreateStatusBar();
  this->ResetActions();
  message_timeout_ = 0;
}


MainWindow::~MainWindow() {
  delete multisnake_;
  delete viewer_;
}

void MainWindow::CreateActions() {
  this->CreateFileMenuActions();
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

  save_as_isotropic_image_ = new QAction(tr("Save as Isotropic Image"), this);
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
  connect(compare_snakes_, SIGNAL(triggered()), this, SLOT(CompareSnakes()));

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

  deform_snakes_in_action_ = new QAction(tr("Deform Snakes in Action"), this);
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
  connect(load_viewpoint_, SIGNAL(triggered()), this, SLOT(LoadViewpoint()));

  save_viewpoint_ = new QAction(tr("Save Viewpoint"), this);
  connect(save_viewpoint_, SIGNAL(triggered()), this, SLOT(SaveViewpoint()));

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

  view_ = menuBar()->addMenu(tr("&View"));
  view_->addAction(toggle_planes_);
  view_->addAction(toggle_mip_);
  view_->addAction(toggle_orientation_marker_);
  view_->addAction(toggle_corner_text_);
  view_->addAction(toggle_bounding_box_);
  view_->addAction(toggle_cube_axes_);
  view_->addAction(toggle_snakes_);
  view_->addAction(toggle_junctions_);

  process_ = menuBar()->addMenu(tr("&Process"));
  process_->addAction(initialize_snakes_);
  process_->addAction(deform_snakes_);
  process_->addAction(deform_snakes_in_action_);
  process_->addAction(deform_one_snake_);
  process_->addAction(cut_snakes_);
  process_->addAction(group_snakes_);

  analysis_ = menuBar()->addMenu(tr("Analysis"));
  actin_cable_submenu_ = analysis_->addMenu(tr("Actin Cable"));
  contractile_ring_submenu_ = analysis_->addMenu(tr("Contractile Ring"));
  fibrin_submenu_ = analysis_->addMenu(tr("Fibrin Network"));
  fibrin_submenu_->addAction(compute_spherical_orientation_);
  droplet_submenu_ = analysis_->addMenu(tr("Droplet"));
  droplet_submenu_->addAction(compute_radial_orientation_);
  droplet_submenu_->addAction(compute_point_density_);
  droplet_submenu_->addAction(compute_curvature_);
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
  toolbar_->addAction(toggle_snakes_);
  toolbar_->addAction(toggle_junctions_);

  toolbar_->addSeparator();
  toolbar_->addAction(initialize_snakes_);
  toolbar_->addAction(deform_snakes_);
  toolbar_->addAction(cut_snakes_);
  toolbar_->addAction(group_snakes_);
  toolbar_->addAction(deform_one_snake_);

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

  toggle_planes_->setEnabled(false);
  toggle_mip_->setEnabled(false);
  toggle_orientation_marker_->setEnabled(false);
  toggle_corner_text_->setEnabled(false);
  toggle_bounding_box_->setEnabled(false);
  toggle_cube_axes_->setEnabled(false);
  toggle_snakes_->setEnabled(false);
  toggle_junctions_->setEnabled(false);

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
  show_analysis_options_->setEnabled(false);

  load_viewpoint_->setEnabled(false);
  save_viewpoint_->setEnabled(false);
  save_snapshot_->setEnabled(false);
}

QString MainWindow::GetLastDirectory(const std::string &filename) {
  QString dir = "..";
  if (!filename.empty()) {
    std::string::size_type pos = filename.find_last_of("/\\");
    dir = QString(filename.substr(0, pos).c_str());
  }
  return dir;
}

void MainWindow::OpenImage() {
  QString dir = this->GetLastDirectory(image_filename_);
  image_filename_ = QFileDialog::getOpenFileName(
      this, tr("Open an image"), dir,
      tr("Image Files (*.tif *.tiff *.mhd *.mha)")).toStdString();
  if (image_filename_.empty()) return;

  this->setWindowTitle(QString("SOAX - ") + image_filename_.c_str());
  multisnake_->LoadImage(image_filename_);
  viewer_->SetupImage(multisnake_->image());
  toggle_planes_->setChecked(true);
  toggle_mip_->setChecked(false);
  toggle_orientation_marker_->setChecked(true);
  toggle_corner_text_->setChecked(true);
  toggle_bounding_box_->setChecked(false);
  toggle_cube_axes_->setChecked(false);
  statusBar()->showMessage(tr("Image loaded."), message_timeout_);

  this->SetParameters();
  analysis_options_dialog_->SetImageCenter(multisnake_->GetImageCenter());

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
  initialize_snakes_->setEnabled(true);
}

void MainWindow::SaveAsIsotropicImage() {
  bool ok;
  double z_spacing = QInputDialog::getDouble(this,
                                             tr("Set inter-slice spacing"),
                                             tr("Z Spacing"),
                                             1.0, 0.1, 10, 4, &ok);
  if (!ok) return;

  if (std::fabs(z_spacing - 1) > kEpsilon) {
    QString dir = this->GetLastDirectory(image_filename_);
    image_filename_ = QFileDialog::getSaveFileName(
        this, tr("Save as Isotropic Image"), dir).toStdString();
    if (image_filename_.empty()) return;
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
  parameter_filename_ = QFileDialog::getOpenFileName(
      this, tr("Open a parameter file"), dir,
      tr("Text Files (*.txt)")).toStdString();
  if (parameter_filename_.empty()) return;
  multisnake_->LoadParameters(parameter_filename_);
  statusBar()->showMessage(tr("Parameters loaded."), message_timeout_);
}

void MainWindow::SaveParameters() {
  QString dir = this->GetLastDirectory(parameter_filename_);
  parameter_filename_ = QFileDialog::getSaveFileName(
      this, tr("Save current parameters"), dir,
      tr("Text Files (*.txt)")).toStdString();
  if (parameter_filename_.empty()) return;
  multisnake_->SaveParameters(parameter_filename_);
  statusBar()->showMessage("Parameters saved.", message_timeout_);
}

void MainWindow::LoadSnakes() {
  QString dir = this->GetLastDirectory(snake_filename_);
  snake_filename_ = QFileDialog::getOpenFileName(
      this, tr("Load Snakes"), dir, "Text (*.txt)").toStdString();
  if (snake_filename_.empty()) return;

  multisnake_->LoadConvergedSnakes(snake_filename_);
  unsigned num_snakes = multisnake_->GetNumberOfConvergedSnakes();
  QString msg = "Number of snakes loaded: " + QString::number(num_snakes);
  statusBar()->showMessage(msg, message_timeout_);
  // multisnake_->PrintSnakes(multisnake_->converged_snakes());
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
  initialize_snakes_->setEnabled(false);
  deform_one_snake_->setEnabled(true);
  toggle_snakes_->setEnabled(true);
  toggle_junctions_->setEnabled(true);
  compute_spherical_orientation_->setEnabled(true);
  compute_radial_orientation_->setEnabled(true);
  compute_point_density_->setEnabled(true);
  compute_curvature_->setEnabled(true);
  show_analysis_options_->setEnabled(true);
}

void MainWindow::SaveSnakes() {
  QString dir = this->GetLastDirectory(snake_filename_);
  snake_filename_ = QFileDialog::getSaveFileName(
      this, tr("Save Snakes"), dir, "Text (*.txt)").toStdString();
  if (snake_filename_.empty()) return;
  multisnake_->SaveSnakes(multisnake_->converged_snakes(), snake_filename_);
  statusBar()->showMessage(tr("Snakes are saved."), message_timeout_);
}

void MainWindow::LoadJFilamentSnakes() {
  QString dir = this->GetLastDirectory(snake_filename_);
  snake_filename_ = QFileDialog::getOpenFileName(
      this, tr("Load JFilament Snakes"), dir, "Text (*.txt)").toStdString();
  if (snake_filename_.empty()) return;
  multisnake_->LoadGroundTruthSnakes(snake_filename_);
  unsigned num_snakes = multisnake_->GetNumberOfComparingSnakes1();
  QString msg = "Number of JFilament snakes loaded: " +
      QString::number(num_snakes);
  statusBar()->showMessage(msg, message_timeout_);
  viewer_->RemoveSnakes();
  viewer_->SetupSnakes(multisnake_->converged_snakes());
  // multisnake_->PrintSnakes(multisnake_->comparing_snakes1());
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
  toggle_snakes_->setEnabled(true);
}

void MainWindow::SaveJFilamentSnakes() {
  QString dir = this->GetLastDirectory(snake_filename_);
  snake_filename_ = QFileDialog::getSaveFileName(
      this, tr("Save JFilament Snakes"), dir, "Text (*.txt)").toStdString();
  if (snake_filename_.empty()) return;
  multisnake_->SaveConvergedSnakesAsJFilamentFormat(snake_filename_);
  statusBar()->showMessage(tr("Snakes are saved in JFilament format"), message_timeout_);
}

void MainWindow::CompareSnakes() {
  QString dir = this->GetLastDirectory(snake_filename_);
  snake_filename_ = QFileDialog::getOpenFileName(
      this, tr("Compare Snakes"), dir, "Text (*.txt)").toStdString();
  if (snake_filename_.empty()) return;
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
  toggle_snakes_->setEnabled(true);
}

void MainWindow::CompareAnotherSnakes() {
  QString dir = this->GetLastDirectory(snake_filename_);
  snake_filename_ = QFileDialog::getOpenFileName(
      this, tr("Compare other snakes"), dir, "Text (*.txt)").toStdString();
  if (snake_filename_.empty()) return;
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
  viewer_->Reset();
  viewer_->Render();
  multisnake_->Reset();
  this->ResetActions();
  open_image_->setEnabled(true);
  setWindowTitle("SOAX");
}

void MainWindow::InitializeSnakes() {
  std::cout << "============ Current Parameters ============" << std::endl;
  multisnake_->WriteParameters(std::cout);
  std::cout << "============================================" << std::endl;

  multisnake_->ScaleImageIntensity();
  multisnake_->ComputeImageGradient();
  multisnake_->InitializeSnakes();
  QString msg = QString::number(multisnake_->GetNumberOfInitialSnakes()) +
      " snakes initialized.";
  statusBar()->showMessage(msg, message_timeout_);

  viewer_->SetupSnakes(multisnake_->initial_snakes());
  toggle_snakes_->setChecked(true);

  initialize_snakes_->setEnabled(false);
  save_snakes_->setEnabled(true);
  toggle_snakes_->setEnabled(true);
  deform_snakes_->setEnabled(true);
  deform_snakes_in_action_->setEnabled(true);
  deform_one_snake_->setEnabled(true);
}

void MainWindow::DeformSnakes() {
  progress_bar_->setMaximum(multisnake_->GetNumberOfInitialSnakes());
  multisnake_->DeformSnakes(progress_bar_);
  unsigned num_snakes = multisnake_->GetNumberOfConvergedSnakes();
  QString msg = "Number of converged snakes: " + QString::number(num_snakes);
  statusBar()->showMessage(msg, message_timeout_);
  viewer_->RemoveSnakes();
  viewer_->SetupSnakes(multisnake_->converged_snakes());
  toggle_snakes_->setChecked(true);
  viewer_->Render();

  deform_snakes_->setEnabled(false);
  toggle_snakes_->setEnabled(true);
  save_snakes_->setEnabled(true);
  save_jfilament_snakes_->setEnabled(true);
  compare_snakes_->setEnabled(true);
  cut_snakes_->setEnabled(true);
  show_analysis_options_->setEnabled(true);
}

void MainWindow::DeformSnakesInAction() {}

void MainWindow::DeformOneSnake() {}

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
}

void MainWindow::ComputeSphericalOrientation() {
}

void MainWindow::ComputeRadialOrientation() {
  QString filename = QFileDialog::getSaveFileName(
      this, tr("Save radial orientation file"), "..", tr("Text (*.txt)"));
  if (filename.isEmpty()) return;

  PointType center;
  analysis_options_dialog_->GetImageCenter(center);
  std::cout << "Image center: " << center << std::endl;
  multisnake_->ComputeRadialOrientation(center, filename.toStdString());
  statusBar()->showMessage(tr("Radial orientation file saved."));
}

void MainWindow::ComputePointDensity() {
  QString filename = QFileDialog::getSaveFileName(
      this, tr("Save snake point density file"), "..", tr("Text (*.txt)"));
  if (filename.isEmpty()) return;
  PointType center;
  analysis_options_dialog_->GetImageCenter(center);
  std::cout << "Image center: " << center << std::endl;
  double radius = analysis_options_dialog_->GetRadius();
  std::cout << "Radius: " << radius << std::endl;
  multisnake_->ComputePointDensity(center, radius, filename.toStdString());
  statusBar()->showMessage(tr("Snake point density file saved."));
}

void MainWindow::ComputeCurvature() {}

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
  multisnake_->set_sigma(parameters_dialog_->GetSigma());
  multisnake_->set_ridge_threshold(parameters_dialog_->GetRidgeThreshold());
  multisnake_->set_foreground(parameters_dialog_->GetForeground());
  multisnake_->set_background(parameters_dialog_->GetBackground());
  Snake::set_background(multisnake_->background() *
                        multisnake_->intensity_scaling());
  Snake::set_desired_spacing(parameters_dialog_->GetSpacing());
  multisnake_->set_initialize_z(parameters_dialog_->InitializeZ());
  Snake::set_minimum_length(parameters_dialog_->GetMinSnakeLength());
  Snake::set_max_iterations(parameters_dialog_->GetMaxIterations());
  Snake::set_change_threshold(parameters_dialog_->GetChangeThreshold());
  Snake::set_check_period(parameters_dialog_->GetCheckPeriod());
  Snake::set_iterations_per_press(parameters_dialog_->GetIterationsPerPress());
  multisnake_->solver_bank()->set_alpha(parameters_dialog_->GetAlpha());
  multisnake_->solver_bank()->set_beta(parameters_dialog_->GetBeta());
  multisnake_->solver_bank()->set_gamma(parameters_dialog_->GetGamma());
  Snake::set_gamma(multisnake_->solver_bank()->gamma());
  Snake::set_external_factor(parameters_dialog_->GetExternalFactor());
  Snake::set_stretch_factor(parameters_dialog_->GetStretchFactor());
  Snake::set_number_of_sectors(parameters_dialog_->GetNumberOfSectors());
  Snake::set_radial_near(parameters_dialog_->GetRadialNear());
  Snake::set_radial_far(parameters_dialog_->GetRadialFar());
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
  viewpoint_filename_ = QFileDialog::getOpenFileName(
      this, tr("Load viewpoint"), dir,
      tr("Camera files (*.cam)")).toStdString();
  if (viewpoint_filename_.empty()) return;
  viewer_->LoadViewpoint(viewpoint_filename_);
  statusBar()->showMessage(tr("Viewpoint loaded."), message_timeout_);
}

void MainWindow::SaveViewpoint() {
  QString dir = this->GetLastDirectory(viewpoint_filename_);
  viewpoint_filename_ = QFileDialog::getSaveFileName(
      this, tr("Save viewpoint"), dir,
      tr("Camera files (*.cam)")).toStdString();
  if (viewpoint_filename_.empty()) return;
  viewer_->SaveViewpoint(viewpoint_filename_);
  statusBar()->showMessage("Viewpoint saved.", message_timeout_);
}

void MainWindow::SaveSnapshot() {
  QString dir = this->GetLastDirectory(snapshot_filename_);
  snapshot_filename_ = QFileDialog::getSaveFileName(
      this, tr("Save snapshot"), dir,
      tr("Image files (*.png *.tif)")).toStdString();
  if (snapshot_filename_.empty()) return;

  viewer_->PrintScreenAsPNGImage(snapshot_filename_);
  viewer_->PrintScreenAsTIFFImage(snapshot_filename_);
  // viewer_->PrintScreenAsVectorImage(snapshot_filename_);
  statusBar()->showMessage("Snapshots saved.", message_timeout_);
}

void MainWindow::AboutSOAX() {
  QMessageBox::about(this, tr("About SOAX"),
                     tr("<h3>SOAX 3.5.0</h3>"
                        "<p>Copyright &copy; 2013 Ting Xu, IDEA Lab, "
                        "Lehigh University "
                        "<p>SOAX extracts curvilinear network structure "
                        "in biomedical images."
                        "This work is supported by NIH grant R01GM098430."));
}

} // namespace soax
