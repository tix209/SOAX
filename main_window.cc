#include <QtGui>
#include "QVTKWidget.h"
#include "main_window.h"
#include "multisnake.h"
#include "viewer.h"

namespace soax {

MainWindow::MainWindow() {
  central_widget_ = new QWidget;
  multisnake_ = new Multisnake;
  viewer_ = new Viewer;
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
  viewer_->RemoveJunctions();
  viewer_->RemoveSnakes();
  viewer_->SetupSnakes(multisnake_->converged_snakes());
  viewer_->SetupSnakes(multisnake_->comparing_snakes1(), 1);
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
  viewer_->SetupSnakes(multisnake_->comparing_snakes1(), 1);
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

void MainWindow::SaveJFilamentSnakes() {}

void MainWindow::CompareSnakes() {}

void MainWindow::CompareAnotherSnakes() {}

void MainWindow::CloseSession() {}

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
