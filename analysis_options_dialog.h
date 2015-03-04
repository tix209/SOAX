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

  void SetImageCenter(const PointType &center);
  void GetImageCenter(PointType &center) const;

  int GetCoarseGraining() const;

  double GetCenterX() const;
  double GetCenterY() const;
  double GetCenterZ() const;
  void SetCenterX(double value);
  void SetCenterY(double value);
  void SetCenterZ(double value);

  unsigned GetRadius() const;
  double GetPixelSize() const;

 public slots:  // NOLINT(whitespace/indent)
  void EnableOKButton();
  void DisableOKButton();

 private:
  QGroupBox *CreatePointDensityGroup();
  QGroupBox *CreateCurvatureGroup();

  QLineEdit *coarse_graining_edit_;
  QLineEdit *center_x_edit_;
  QLineEdit *center_y_edit_;
  QLineEdit *center_z_edit_;
  QLineEdit *radius_edit_;
  QLineEdit *pixel_size_edit_;
  QDialogButtonBox *button_box_;

  DISALLOW_COPY_AND_ASSIGN(AnalysisOptionsDialog);
};

}  // namespace soax

#endif  // ANALYSIS_OPTIONS_DIALOG_H_
