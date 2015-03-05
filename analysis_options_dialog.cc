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
  layout->addWidget(this->CreateGeneralGroup());
  layout->addWidget(this->CreateCurvatureGroup());
  layout->addWidget(this->CreatePointDensityGroup());

  button_box_ = new QDialogButtonBox(QDialogButtonBox::Ok |
                                     QDialogButtonBox::Cancel);

  layout->addWidget(button_box_);
  setLayout(layout);
  setWindowTitle(tr("Analysis Options"));
  connect(button_box_, SIGNAL(accepted()), this, SLOT(accept()));
  connect(button_box_, SIGNAL(rejected()), this, SLOT(reject()));
}

double AnalysisOptionsDialog::GetPixelSize() const {
  return pixel_size_edit_->text().toDouble();
}

int AnalysisOptionsDialog::GetCoarseGraining() const {
  return coarse_graining_edit_->text().toInt();
}

void AnalysisOptionsDialog::GetImageCenter(PointType &center) const {
  center[0] = center_x_edit_->text().toDouble();
  center[1] = center_y_edit_->text().toDouble();
  center[2] = center_z_edit_->text().toDouble();
}

void AnalysisOptionsDialog::SetImageCenter(const PointType &center) {
  center_x_edit_->setText(QString::number(center[0]));
  center_y_edit_->setText(QString::number(center[1]));
  center_z_edit_->setText(QString::number(center[2]));
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

void AnalysisOptionsDialog::SetRadius(double r) {
  radius_edit_->setText(QString::number(r));
}

void AnalysisOptionsDialog::DisableOKButton() {
  button_box_->button(QDialogButtonBox::Ok)->setEnabled(false);
}

void AnalysisOptionsDialog::EnableOKButton() {
  button_box_->button(QDialogButtonBox::Ok)->setEnabled(true);
}

QGroupBox * AnalysisOptionsDialog::CreateGeneralGroup() {
  QGroupBox *gb = new QGroupBox(tr("General"));
  pixel_size_edit_ = new QLineEdit("1.0");
  connect(pixel_size_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));
  QLabel *label = new QLabel(tr("Pixel Size (um)"));
  QHBoxLayout *hbox = new QHBoxLayout;
  hbox->addWidget(label);
  hbox->addWidget(pixel_size_edit_);
  hbox->addStretch();
  gb->setLayout(hbox);
  return gb;
}

QGroupBox * AnalysisOptionsDialog::CreateCurvatureGroup() {
  QGroupBox *gb = new QGroupBox(tr("Curvature"));
  coarse_graining_edit_ = new QLineEdit("8");
  connect(coarse_graining_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));
  QLabel *label = new QLabel(tr("Coarse Graining (snake points)"));
  QHBoxLayout *hbox = new QHBoxLayout;
  hbox->addWidget(label);
  hbox->addWidget(coarse_graining_edit_);
  hbox->addStretch();
  gb->setLayout(hbox);
  return gb;
}

QGroupBox * AnalysisOptionsDialog::CreatePointDensityGroup() {
  QGroupBox *gb = new QGroupBox(tr("Spherical Confinement"));
  center_x_edit_ = new QLineEdit("0");
  center_y_edit_ = new QLineEdit("0");
  center_z_edit_ = new QLineEdit("0");
  radius_edit_ = new QLineEdit("0");

  connect(center_x_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(center_y_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(center_z_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(radius_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));

  QLabel *label_x = new QLabel(tr("Center X (pixels)"));
  QLabel *label_y = new QLabel(tr("Center Y (pixels)"));
  QLabel *label_z = new QLabel(tr("Center Z (pixels)"));
  QLabel *label_r = new QLabel(tr("Radius (pixels)"));
  QHBoxLayout *hbox1 = new QHBoxLayout;
  hbox1->addWidget(label_x);
  hbox1->addWidget(center_x_edit_);
  hbox1->addStretch();
  QHBoxLayout *hbox2 = new QHBoxLayout;
  hbox2->addWidget(label_y);
  hbox2->addWidget(center_y_edit_);
  hbox2->addStretch();
  QHBoxLayout *hbox3 = new QHBoxLayout;
  hbox3->addWidget(label_z);
  hbox3->addWidget(center_z_edit_);
  hbox3->addStretch();
  QHBoxLayout *hbox4 = new QHBoxLayout;
  hbox4->addWidget(label_r);
  hbox4->addWidget(radius_edit_);
  hbox4->addStretch();
  QVBoxLayout *vbox = new QVBoxLayout;
  vbox->addLayout(hbox1);
  vbox->addLayout(hbox2);
  vbox->addLayout(hbox3);
  vbox->addLayout(hbox4);
  gb->setLayout(vbox);
  return gb;
}


}  // namespace soax
