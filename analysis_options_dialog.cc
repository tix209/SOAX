/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file implements the analysis dialog for SOAX.
 */

#include "./analysis_options_dialog.h"
#include <QtWidgets>  // NOLINT(build/include_order)

namespace soax {
AnalysisOptionsDialog::AnalysisOptionsDialog(QWidget *parent) :
    QDialog(parent) {
  QVBoxLayout *layout = new QVBoxLayout;
  layout->addWidget(this->CreateGeneralGroup());
  layout->addWidget(this->CreateCurvatureGroup());
  layout->addWidget(this->CreateSphericalConfinementGroup());

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

bool AnalysisOptionsDialog::GetInsideRatio(double *ratio) const {
  bool ok = false;
  *ratio = inside_ratio_edit_->text().toDouble(&ok);
  return ok;
}

bool AnalysisOptionsDialog::ExcludeBoundaryChecked() const {
  return exclude_boundary_check_->isChecked();
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
  exclude_boundary_check_ = new QCheckBox(tr(""));
  connect(pixel_size_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(exclude_boundary_check_, SIGNAL(stateChanged(int)),
          this, SLOT(EnableOKButton()));

  QHBoxLayout *hbox = new QHBoxLayout;
  hbox->addWidget(exclude_boundary_check_);
  hbox->addStretch();

  QFormLayout *form = new QFormLayout;
  form->addRow(tr("Pixel Size (um)"), pixel_size_edit_);
  form->addRow(tr("Exclude points near image boundary"), hbox);

  gb->setLayout(form);
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

QGroupBox * AnalysisOptionsDialog::CreateSphericalConfinementGroup() {
  QGroupBox *gb = new QGroupBox(tr("Spherical Confinement"));
  for (int i = 0; i < kDimension; i++) {
    center_edit_[i] = new QLineEdit("0.0");
    connect(center_edit_[i], SIGNAL(textChanged(const QString &)),
            this, SLOT(EnableOKButton()));
  }
  radius_edit_ = new QLineEdit("0");
  inside_ratio_edit_ = new QLineEdit("1.0");

  connect(radius_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));

  QFormLayout *form  = new QFormLayout;
  form->addRow(tr("Center X (pixels)"), center_edit_[0]);
  form->addRow(tr("Center Y (pixels)"), center_edit_[1]);
  form->addRow(tr("Center Z (pixels)"), center_edit_[2]);
  form->addRow(tr("Radius (pixels)"), radius_edit_);
  form->addRow(tr("Inside Ratio"), inside_ratio_edit_);

  gb->setLayout(form);
  return gb;
}


}  // namespace soax
