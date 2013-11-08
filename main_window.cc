#include <QtGui>
// #include "QVTKWidget.h"
#include "main_window.h"

namespace soax {

MainWindow::MainWindow() {
  central_widget_ = new QWidget;

  QHBoxLayout* layout = new QHBoxLayout;
  // layout->addWidget(viewer_->qvtk());
  central_widget_->setLayout(layout);
  setCentralWidget(central_widget_);
  setWindowIcon(QIcon(":/icon/letter-x.png"));

  this->CreateActions();
  this->CreateMenus();

}


MainWindow::~MainWindow() {
  // delete viewer_;
  // delete solver_bank_;
}

void MainWindow::CreateActions() {
  this->CreateHelpMenuActions();
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
  help_menu_ = menuBar()->addMenu(tr("&Help"));
  help_menu_->addAction(about_soax_action_);
  help_menu_->addAction(about_qt_action_);
}


void MainWindow::AboutSOAX() {
  QMessageBox::about(this, tr("About SOAX"),
                     tr("<h3>SOAX 3.5.0</h3>"
                        "<p>Copyright &copy; 2013 Ting Xu, IDEA Lab, "
                        "Lehigh University "
                        "<p>SOAX is an application that "
                        "extracts curvilinear network structure "
                        "from 3D biomedical images. "
                        "This work is supported by NIH, grant R01GM098430."));
}

} // namespace soax
