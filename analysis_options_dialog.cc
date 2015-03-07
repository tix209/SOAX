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

bool AnalysisOptionsDialog::GetPixelSize(double *pixel_size) const {
  bool ok = false;
  *pixel_size = pixel_size_edit_->text().toDouble(&ok);
  return ok;
}

bool AnalysisOptionsDialog::GetCoarseGraining(int *coarse_graining) const {
  bool ok = false;
  *coarse_graining = coarse_graining_edit_->text().toInt(&ok);
  return ok;
}

bool AnalysisOptionsDialog::GetImageCenter(PointType *center) const {
  bool ok = false;
  for (int i = 0; i < kDimension; i++) {
    (*center)[i] = center_edit_[i]->text().toDouble(&ok);
    if (!ok) return false;
  }
  return true;
}

void AnalysisOptionsDialog::SetImageCenter(const PointType &center) {
  for (int i = 0; i < kDimension; i++) {
    center_edit_[i]->setText(QString::number(center[i]));
  }
}

bool AnalysisOptionsDialog::GetRadius(double *radius) const {
  bool ok = false;
  *radius = radius_edit_->text().toDouble(&ok);
  return ok;
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
  for (int i = 0; i < kDimension; i++) {
    center_edit_[i] = new QLineEdit("0.0");
    connect(center_edit_[i], SIGNAL(textChanged(const QString &)),
            this, SLOT(EnableOKButton()));
  }
  radius_edit_ = new QLineEdit("0");

  connect(radius_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));

  QLabel *label_x = new QLabel(tr("Center X (pixels)"));
  QLabel *label_y = new QLabel(tr("Center Y (pixels)"));
  QLabel *label_z = new QLabel(tr("Center Z (pixels)"));
  QLabel *label_r = new QLabel(tr("Radius (pixels)"));
  QHBoxLayout *hbox1 = new QHBoxLayout;
  hbox1->addWidget(label_x);
  hbox1->addWidget(center_edit_[0]);
  hbox1->addStretch();
  QHBoxLayout *hbox2 = new QHBoxLayout;
  hbox2->addWidget(label_y);
  hbox2->addWidget(center_edit_[1]);
  hbox2->addStretch();
  QHBoxLayout *hbox3 = new QHBoxLayout;
  hbox3->addWidget(label_z);
  hbox3->addWidget(center_edit_[2]);
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
