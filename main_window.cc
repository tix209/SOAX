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
  open_image_action_ = new QAction(tr("&Open Image"), this);
  open_image_action_->setShortcut(QKeySequence::Open);
  open_image_action_->setIcon(QIcon(":/icon/Open.png"));
  connect(open_image_action_, SIGNAL(triggered()),
          this, SLOT(OpenImage()));
}

void MainWindow::CreateViewMenuActions() {
  toggle_planes_action_ = new QAction(tr("Slice Planes"), this);
  toggle_planes_action_->setIcon(QIcon(":/icon/Picture.png"));
  toggle_planes_action_->setCheckable(true);
  connect(toggle_planes_action_, SIGNAL(toggled(bool)),
          viewer_, SLOT(ToggleSlicePlanes(bool)));
}

void MainWindow::CreateHelpMenuActions() {
  about_soax_action_ = new QAction(tr("About SOAX"), this);
  about_soax_action_->setStatusTip(tr("About SOAX"));
  connect(about_soax_action_, SIGNAL(triggered()), this, SLOT(AboutSOAX()));

  about_qt_action_ = new QAction(tr("About Qt"), this);
  about_qt_action_->setStatusTip(tr("About Qt"));
  connect(about_qt_action_, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void MainWindow::CreateMenus() {
  file_menu_ = menuBar()->addMenu(tr("&File"));
  file_menu_->addAction(open_image_action_);

  view_menu_ = menuBar()->addMenu(tr("&View"));
  view_menu_->addAction(toggle_planes_action_);

  help_menu_ = menuBar()->addMenu(tr("&Help"));
  help_menu_->addAction(about_soax_action_);
  help_menu_->addAction(about_qt_action_);
}

void MainWindow::CreateToolBar() {
  toolbar_ = addToolBar(tr("shortcuts"));
  toolbar_->addAction(open_image_action_);

  toolbar_->addSeparator();
  toolbar_->addAction(toggle_planes_action_);
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
  this->setWindowTitle(QString("soax3D - ") + image_filename_.c_str());
  multisnake_->LoadImage(image_filename_);
  viewer_->SetupImage(multisnake_->image());
  // viewer_->DisplayOrientationMarker();
  // viewer_->DisplayBoundingBox();
  // viewer_->DisplayUpperLeftCornerInfo();
  toggle_planes_action_->setChecked(true);
  // toggle_volume_action_->setChecked(false);
  viewer_->Render();

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
