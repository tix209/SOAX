/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file defines the analysis dialog for SOAX.
 */

#ifndef ANALYSIS_OPTIONS_DIALOG_H_
#define ANALYSIS_OPTIONS_DIALOG_H_

#include <QDialog>
#include "./global.h"
class QLineEdit;
class QGroupBox;
class QDialogButtonBox;

namespace soax {

class AnalysisOptionsDialog : public QDialog {
  Q_OBJECT

 public:
  explicit AnalysisOptionsDialog(QWidget * parent = NULL);

  bool GetPixelSize(double *pixel_size) const;
  bool GetCoarseGraining(int *coarse_graining) const;

  bool GetImageCenter(PointType *center) const;
  void SetImageCenter(const PointType &center);

  bool GetRadius(double *radius) const;
  void SetRadius(double r);

 public slots:  // NOLINT(whitespace/indent)
  void EnableOKButton();
  void DisableOKButton();

 private:
  QGroupBox *CreateGeneralGroup();
  QGroupBox *CreateCurvatureGroup();
  QGroupBox *CreatePointDensityGroup();

  QLineEdit *coarse_graining_edit_;
  QLineEdit *center_edit_[kDimension];
  QLineEdit *radius_edit_;
  QLineEdit *pixel_size_edit_;
  QDialogButtonBox *button_box_;

  DISALLOW_COPY_AND_ASSIGN(AnalysisOptionsDialog);
};

}  // namespace soax

#endif  // ANALYSIS_OPTIONS_DIALOG_H_
