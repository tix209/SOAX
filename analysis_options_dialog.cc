/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file implements the analysis dialog for SOAX.
 */

#include "./analysis_options_dialog.h"
#include <QtGui>  // NOLINT(build/include_order)

namespace soax {
AnalysisOptionsDialog::AnalysisOptionsDialog(QWidget *parent) :
    QDialog(parent) {
  QVBoxLayout *layout = new QVBoxLayout;
  layout->addWidget(this->CreatePointDensityGroup());
  layout->addWidget(this->CreateCurvatureGroup());

  button_box_ = new QDialogButtonBox(QDialogButtonBox::Ok |
                                     QDialogButtonBox::Cancel);

  layout->addWidget(button_box_);
  setLayout(layout);
  setWindowTitle(tr("Analysis Options"));
  connect(button_box_, SIGNAL(accepted()), this, SLOT(accept()));
  connect(button_box_, SIGNAL(rejected()), this, SLOT(reject()));
}

void AnalysisOptionsDialog::SetImageCenter(const PointType &center) {
  center_x_edit_->setText(QString::number(center[0]));
  center_y_edit_->setText(QString::number(center[1]));
  center_z_edit_->setText(QString::number(center[2]));
}

void AnalysisOptionsDialog::GetImageCenter(PointType &center) const {
  center[0] = center_x_edit_->text().toDouble();
  center[1] = center_y_edit_->text().toDouble();
  center[2] = center_z_edit_->text().toDouble();
}

QGroupBox * AnalysisOptionsDialog::CreateCurvatureGroup() {
  QGroupBox *gb = new QGroupBox(tr("Curvature (1/pixel)"));
  coarse_graining_edit_ = new QLineEdit("8");

  connect(coarse_graining_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));

  QLabel *label = new QLabel(tr("Coarse graining (SOAC points)"));

  QHBoxLayout *layout = new QHBoxLayout;
  layout->addWidget(label);
  layout->addWidget(coarse_graining_edit_);
  layout->addStretch();
  gb->setLayout(layout);
  return gb;
}

QGroupBox * AnalysisOptionsDialog::CreatePointDensityGroup() {
  QGroupBox *gb = new QGroupBox(
      tr("Spherical, Radial Orientation && Point Density"));
  center_x_edit_ = new QLineEdit("0");
  center_y_edit_ = new QLineEdit("0");
  center_z_edit_ = new QLineEdit("0");
  radius_edit_ = new QLineEdit("100");
  pixel_size_edit_ = new QLineEdit("1.0");

  connect(center_x_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(center_y_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(center_z_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(radius_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(pixel_size_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));

  QLabel *center_label = new QLabel(tr("Center x, y, z (pixels)"));
  QLabel *radius_label = new QLabel(tr("Radius (pixels) "));
  QLabel *pixel_size_label = new QLabel(tr("Pixel size (um)"));

  QHBoxLayout *hlayout1 = new QHBoxLayout;
  hlayout1->addWidget(center_label);
  hlayout1->addWidget(center_x_edit_);
  hlayout1->addWidget(center_y_edit_);
  hlayout1->addWidget(center_z_edit_);
  hlayout1->addStretch();

  QHBoxLayout *hlayout4 = new QHBoxLayout;
  hlayout4->addWidget(radius_label);
  hlayout4->addWidget(radius_edit_);
  hlayout4->addWidget(pixel_size_label);
  hlayout4->addWidget(pixel_size_edit_);
  hlayout4->addStretch();

  QVBoxLayout *vlayout  = new QVBoxLayout;
  vlayout->addLayout(hlayout1);
  vlayout->addLayout(hlayout4);

  gb->setLayout(vlayout);
  return gb;
}

void AnalysisOptionsDialog::DisableOKButton() {
  button_box_->button(QDialogButtonBox::Ok)->setEnabled(false);
}

void AnalysisOptionsDialog::EnableOKButton() {
  button_box_->button(QDialogButtonBox::Ok)->setEnabled(true);
}

int AnalysisOptionsDialog::GetCoarseGraining() const {
  return coarse_graining_edit_->text().toInt();
}

double AnalysisOptionsDialog::GetCenterX() const {
  return center_x_edit_->text().toDouble();
}

double AnalysisOptionsDialog::GetCenterY() const {
  return center_y_edit_->text().toDouble();
}

double AnalysisOptionsDialog::GetCenterZ() const {
  return center_z_edit_->text().toDouble();
}

void AnalysisOptionsDialog::SetCenterX(double value) {
  QString s;
  s.setNum(value);
  center_x_edit_->setText(s);
}

void AnalysisOptionsDialog::SetCenterY(double value) {
  QString s;
  s.setNum(value);
  center_y_edit_->setText(s);
}

void AnalysisOptionsDialog::SetCenterZ(double value) {
  QString s;
  s.setNum(value);
  center_z_edit_->setText(s);
}

unsigned AnalysisOptionsDialog::GetRadius() const {
  return radius_edit_->text().toUInt();
}

double AnalysisOptionsDialog::GetPixelSize() const {
  return pixel_size_edit_->text().toDouble();
}

}  // namespace soax
