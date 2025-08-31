/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file implements SOAX as a Qt application.
 */

#include <QApplication>
#include "./main_window.h"
#include <vtkRenderingOpenGL2Module.h>
#include <vtkRenderingVolumeOpenGL2Module.h>
#include <vtkAutoInit.h>
#include "QVTKOpenGLWidget.h"
#include "vtkGenericOpenGLRenderWindow.h"
#include <QSurfaceFormat>

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2);

int main(int argc, char **argv) {
#ifndef __APPLE__
  vtkOpenGLRenderWindow::SetGlobalMaximumNumberOfMultiSamples(0);
#endif

  QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());
  QApplication app(argc, argv);

#ifdef __APPLE__
  app.setAttribute(Qt::AA_UseHighDpiPixmaps);
#endif

  soax::MainWindow window;
  window.resize(1024, 768);
  window.show();
  return app.exec();
}
