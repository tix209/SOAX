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
}


MainWindow::~MainWindow() {
  delete multisnake_;
  delete viewer_;
}

void MainWindow::CreateActions() {
  this->CreateFileMenuActions();
  this->CreateViewMenuActions();
  this->CreateHelpMenuActions();
}

void MainWindow::CreateFileMenuActions() {
  open_image_ = new QAction(tr("&Open Image"), this);
  open_image_->setShortcut(QKeySequence::Open);
  open_image_->setIcon(QIcon(":/icon/Open.png"));
  connect(open_image_, SIGNAL(triggered()),
          this, SLOT(OpenImage()));
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

  view_ = menuBar()->addMenu(tr("&View"));
  view_->addAction(toggle_planes_);
  view_->addAction(toggle_mip_);
  view_->addAction(toggle_orientation_marker_);
  view_->addAction(toggle_corner_text_);
  view_->addAction(toggle_bounding_box_);
  view_->addAction(toggle_cube_axes_);

  help_ = menuBar()->addMenu(tr("&Help"));
  help_->addAction(about_soax_);
  help_->addAction(about_qt_);
}

void MainWindow::CreateToolBar() {
  toolbar_ = addToolBar(tr("shortcuts"));
  toolbar_->addAction(open_image_);

  toolbar_->addSeparator();
  toolbar_->addAction(toggle_planes_);
  toolbar_->addAction(toggle_mip_);
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
}

void MainWindow::AboutSOAX() {
  QMessageBox::about(this, tr("About SOAX"),
                     tr("<h3>SOAX 3.5.0</h3>"
                        "<p>Copyright &copy; 2013 Ting Xu, IDEA Lab, "
                        "Lehigh University "
                        "<p>SOAX extracts curvilinear network structure "
                        "in biomedical images."
                        "This work is supported by NIH, grant R01GM098430."));
}

} // namespace soax
