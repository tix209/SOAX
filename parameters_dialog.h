/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file defines the parameter dialog for SOAX.
 */


#ifndef PARAMETERS_DIALOG_H_
#define PARAMETERS_DIALOG_H_


#include <QDialog>
#include <QDialogButtonBox>
#include <QLineEdit>
#include <QCheckBox>
#include "./global.h"
class QGroupBox;

namespace soax {

class Multisnake;


class ParametersDialog : public QDialog {
  Q_OBJECT

 public:
  explicit ParametersDialog(QWidget *parent = NULL);

  double GetIntensityScaling() {
    return intensity_scaling_edit_->text().toDouble();
  }
  double GetSigma() {return sigma_edit_->text().toDouble();}
  double GetRidgeThreshold() {
    return ridge_threshold_edit_->text().toDouble();
  }
  unsigned GetForeground() {return foreground_edit_->text().toUInt();}
  unsigned GetBackground() {return background_edit_->text().toUInt();}
  double GetSpacing() {return spacing_edit_->text().toDouble();}
  bool InitializeZ() {return initialize_z_check_->isChecked();}
  unsigned GetMinSnakeLength() {
    return min_snake_length_edit_->text().toDouble();
  }
  unsigned GetMaxIterations() {
    return max_iterations_edit_->text().toUInt();
  }
  double GetChangeThreshold() {
    return change_threshold_edit_->text().toDouble();
  }
  unsigned GetCheckPeriod() {return check_period_edit_->text().toUInt();}
  unsigned GetIterationsPerPress() {
    return iterations_per_press_edit_->text().toUInt();
  }
  double GetAlpha() {return alpha_edit_->text().toDouble();}
  double GetBeta() {return beta_edit_->text().toDouble();}
  double GetGamma() {return gamma_edit_->text().toDouble();}
  double GetExternalFactor() {
    return external_factor_edit_->text().toDouble();
  }
  double GetStretchFactor() {return stretch_factor_edit_->text().toDouble();}
  unsigned GetNumberOfSectors() {
    return number_of_sectors_edit_->text().toInt();
  }
  unsigned GetRadialNear() {return radial_near_edit_->text().toInt();}
  unsigned GetRadialFar() {return radial_far_edit_->text().toInt();}
  double GetZSpacing() {return z_spacing_edit_->text().toDouble();}
  unsigned GetDelta() {return delta_edit_->text().toUInt();}
  double GetOverlapThreshold() {
    return overlap_threshold_edit_->text().toDouble();
  }
  double GetGroupingDistanceThreshold() {
    return grouping_distance_threshold_edit_->text().toDouble();
  }
  unsigned GetGroupingDelta() {return grouping_delta_edit_->text().toUInt();}
  double GetDirectionThreshold() {
    return direction_threshold_edit_->text().toDouble();
  }
  bool DampZ() {return damp_z_check_->isChecked();}

  void SetCurrentParameters(Multisnake *ms);

 public slots:  // NOLINT(whitespace/indent)
  void EnableOKButton();
  void DisableOKButton();

 private:
  QGroupBox * CreateSnakeGroup();


  QDialogButtonBox *button_box_;
  QLineEdit *intensity_scaling_edit_;
  QLineEdit *sigma_edit_;
  QLineEdit *ridge_threshold_edit_;
  QLineEdit *foreground_edit_;
  QLineEdit *background_edit_;
  QLineEdit *spacing_edit_;
  QLineEdit *min_snake_length_edit_;
  QLineEdit *max_iterations_edit_;
  QLineEdit *change_threshold_edit_;
  QLineEdit *check_period_edit_;
  QLineEdit *iterations_per_press_edit_;

  QLineEdit *alpha_edit_;
  QLineEdit *beta_edit_;
  QLineEdit *gamma_edit_;
  QLineEdit *external_factor_edit_;
  QLineEdit *stretch_factor_edit_;
  QLineEdit *number_of_sectors_edit_;
  QLineEdit *radial_near_edit_;
  QLineEdit *radial_far_edit_;
  QLineEdit *z_spacing_edit_;
  QLineEdit *delta_edit_;
  QLineEdit *overlap_threshold_edit_;
  QLineEdit *grouping_distance_threshold_edit_;
  QLineEdit *grouping_delta_edit_;
  QLineEdit *direction_threshold_edit_;

  QCheckBox *initialize_z_check_;
  QCheckBox *damp_z_check_;

  DISALLOW_COPY_AND_ASSIGN(ParametersDialog);
};

}  // namespace soax

#endif  // PARAMETERS_DIALOG_H_
