/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file implements SOAX as a Qt application.
 */

#include <QApplication>
#include "./main_window.h"

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  soax::MainWindow window;
  window.resize(1024, 768);
  window.show();
  return app.exec();
}
