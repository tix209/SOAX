#ifndef SOAX_MAINWINDOW_H_
#define SOAX_MAINWINDOW_H_

#include <QMainWindow>
#include "global.h"


namespace soax {


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
  // void OpenImage();

  void AboutSOAX();

 private:
  void CreateActions();
  void CreateHelpMenuActions();

  void CreateMenus();


  QWidget *central_widget_;

  // Actions in Help menu
  QAction *about_soax_action_;
  QAction *about_qt_action_;

  // Menus
  QMenu *help_menu_;


  DISALLOW_COPY_AND_ASSIGN(MainWindow);
};

} // namespace soax

#endif // SOAX_MAINWINDOW_H_
