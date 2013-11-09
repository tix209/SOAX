#ifndef SOAX_MAINWINDOW_H_
#define SOAX_MAINWINDOW_H_

#include <QMainWindow>
#include "global.h"


namespace soax {

class Multisnake;
class Viewer;

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

  void AboutSOAX();

 private:
  void CreateActions();
  void CreateFileMenuActions();
  void CreateViewMenuActions();
  void CreateHelpMenuActions();

  void CreateMenus();
  void CreateToolBar();

  QString GetLastDirectory(const std::string &filename);


  QWidget *central_widget_;

  // Actions in File menu
  QAction *open_image_action_;

  // Actions in View menu
  QAction *toggle_planes_action_;

  // Actions in Help menu
  QAction *about_soax_action_;
  QAction *about_qt_action_;

  // Menus
  QMenu *file_menu_;
  QMenu *view_menu_;
  QMenu *help_menu_;

  QToolBar *toolbar_;

  // Complete image file path and name
  std::string image_filename_;
  Multisnake *multisnake_;
  Viewer *viewer_;

  DISALLOW_COPY_AND_ASSIGN(MainWindow);
};

} // namespace soax

#endif // SOAX_MAINWINDOW_H_
