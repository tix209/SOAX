#include <QtGui>
#include "parameters_dialog.h"
#include "multisnake.h"
#include "snake.h"
#include "solver_bank.h"

namespace soax {
ParametersDialog::ParametersDialog(QWidget *parent) : QDialog(parent) {
  QVBoxLayout *layout = new QVBoxLayout;
  //  QGridLayout *layout = new QGridLayout;
  layout->addWidget(this->CreateSnakeGroup());

  button_box_ = new QDialogButtonBox(QDialogButtonBox::Ok |
                                     QDialogButtonBox::Cancel);
  this->DisableOKButton();
  layout->addWidget(button_box_);

  setLayout(layout);
  setWindowTitle(tr("Parameter Settings"));

  connect(button_box_, SIGNAL(accepted()), this, SLOT(accept()));
  connect(button_box_, SIGNAL(rejected()), this, SLOT(reject()));
}

void ParametersDialog::SetCurrentParameters(Multisnake *ms) {
  intensity_scaling_edit_->setText(QString::number(ms->intensity_scaling()));
  sigma_edit_->setText(QString::number(ms->sigma()));
  ridge_threshold_edit_->setText(QString::number(ms->ridge_threshold()));
  foreground_edit_->setText(QString::number(ms->foreground()));
  background_edit_->setText(QString::number(ms->background()));
  spacing_edit_->setText(QString::number(Snake::desired_spacing()));
  min_snake_length_edit_->setText(QString::number(Snake::minimum_length()));
  max_iterations_edit_->setText(QString::number(Snake::max_iterations()));
  change_threshold_edit_->setText(
      QString::number(Snake::change_threshold()));
  check_period_edit_->setText(QString::number(Snake::check_period()));
  iterations_per_press_edit_->setText(
      QString::number(Snake::iterations_per_press()));
  alpha_edit_->setText(QString::number(ms->solver_bank()->alpha()));
  beta_edit_->setText(QString::number(ms->solver_bank()->beta()));
  gamma_edit_->setText(QString::number(ms->solver_bank()->gamma()));
  external_factor_edit_->setText(QString::number(Snake::external_factor()));
  stretch_factor_edit_->setText(QString::number(Snake::stretch_factor()));
  number_of_sectors_edit_->setText(
      QString::number(Snake::number_of_sectors()));
  radial_near_edit_->setText(QString::number(Snake::radial_near()));
  radial_far_edit_->setText(QString::number(Snake::radial_far()));
  delta_edit_->setText(QString::number(Snake::delta()));
  overlap_threshold_edit_->setText(
      QString::number(Snake::overlap_threshold()));
  grouping_distance_threshold_edit_->setText(
      QString::number(Snake::grouping_distance_threshold()));
  grouping_delta_edit_->setText(QString::number(Snake::grouping_delta()));
  direction_threshold_edit_->setText(
      QString::number(Snake::direction_threshold()));
  initialize_z_check_->setChecked(ms->initialize_z());
  damp_z_check_->setChecked(Snake::damp_z());
}

void ParametersDialog::EnableOKButton() {
  button_box_->button(QDialogButtonBox::Ok)->setEnabled(true);
}

void ParametersDialog::DisableOKButton() {
  button_box_->button(QDialogButtonBox::Ok)->setEnabled(false);
}

QGroupBox * ParametersDialog::CreateSnakeGroup() {
  QGroupBox *gp = new QGroupBox("");
  intensity_scaling_edit_ = new QLineEdit("0.0");
  sigma_edit_ = new QLineEdit("0.0");
  ridge_threshold_edit_ = new QLineEdit("0.0");
  foreground_edit_ = new QLineEdit("0");
  background_edit_ = new QLineEdit("0");
  spacing_edit_ = new QLineEdit("1.0");
  min_snake_length_edit_ = new QLineEdit("0.0");
  max_iterations_edit_ = new QLineEdit("0");
  change_threshold_edit_ = new QLineEdit("0.0");
  check_period_edit_ = new QLineEdit("0");
  iterations_per_press_edit_ = new QLineEdit("100");
  alpha_edit_ = new QLineEdit("0.0");
  beta_edit_ = new QLineEdit("0.0");
  gamma_edit_ = new QLineEdit("0.0");
  external_factor_edit_ = new QLineEdit("0.0");
  stretch_factor_edit_ = new QLineEdit("0.0");
  number_of_sectors_edit_ = new QLineEdit("0");
  radial_near_edit_ = new QLineEdit("0");
  radial_far_edit_ = new QLineEdit("0");
  delta_edit_ = new QLineEdit("0");
  overlap_threshold_edit_ = new QLineEdit("0.0");
  grouping_distance_threshold_edit_ = new QLineEdit("0.0");
  grouping_delta_edit_ = new QLineEdit("0");
  direction_threshold_edit_ = new QLineEdit("0.0");

  initialize_z_check_ = new QCheckBox(tr("Init z"));
  initialize_z_check_->setChecked(false);
  damp_z_check_ = new QCheckBox(tr("Damp z"));
  damp_z_check_->setChecked(false);

  connect(intensity_scaling_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(sigma_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(ridge_threshold_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(foreground_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(background_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(spacing_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(min_snake_length_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(max_iterations_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(change_threshold_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(check_period_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(iterations_per_press_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(alpha_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(beta_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(gamma_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(external_factor_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(stretch_factor_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(number_of_sectors_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(radial_near_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(radial_far_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(delta_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(overlap_threshold_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(grouping_distance_threshold_edit_,
          SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(grouping_delta_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(direction_threshold_edit_, SIGNAL(textEdited(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(initialize_z_check_, SIGNAL(stateChanged(int)),
          this, SLOT(EnableOKButton()));
  connect(damp_z_check_, SIGNAL(stateChanged(int)),
          this, SLOT(EnableOKButton()));

  QFormLayout *layout_left  = new QFormLayout;
  layout_left->addRow(tr("Image Intensity Scaling"),
                      intensity_scaling_edit_);
  layout_left->addRow(tr("Gaussian Standard Deviation (pixels)"),
                      sigma_edit_);
  layout_left->addRow(tr("Ridge Threshold"), ridge_threshold_edit_);
  layout_left->addRow(tr("Foreground"), foreground_edit_);
  layout_left->addRow(tr("Background"), background_edit_);
  layout_left->addRow(tr("Snake Point Spacing (pixels) "), spacing_edit_);
  layout_left->addRow(tr("Minimum Snake Length (pixels)"),
                      min_snake_length_edit_);
  layout_left->addRow(tr("Maximum Iterations"), max_iterations_edit_);
  layout_left->addRow(tr("Change Threshold (pixels)"),
                      change_threshold_edit_);
  layout_left->addRow(tr("Check Period"), check_period_edit_);
  layout_left->addRow(tr("Iterations per press"),
                      iterations_per_press_edit_);
  layout_left->addRow(tr(""), initialize_z_check_);
  layout_left->addRow(tr(""), damp_z_check_);

  QFormLayout *layout_right  = new QFormLayout;
  layout_right->addRow(tr("Alpha"), alpha_edit_);
  layout_right->addRow(tr("Beta"), beta_edit_);
  layout_right->addRow(tr("Gamma"), gamma_edit_);
  layout_right->addRow(tr("External Factor"), external_factor_edit_);
  layout_right->addRow(tr("Stretch Factor"), stretch_factor_edit_);
  layout_right->addRow(tr("Number of Sectors"), number_of_sectors_edit_);
  layout_right->addRow(tr("Radial Near (pixels)"), radial_near_edit_);
  layout_right->addRow(tr("Radial Far (pixels)"), radial_far_edit_);
  layout_right->addRow(tr("Delta"), delta_edit_);
  layout_right->addRow(tr("Overlap Threshold (pixels)"),
                       overlap_threshold_edit_);
  layout_right->addRow(tr("Grouping Distance Threshold (pixels)"),
                       grouping_distance_threshold_edit_);
  layout_right->addRow(tr("Grouping Delta"), grouping_delta_edit_);
  layout_right->addRow(tr("Direction Threshold (radians)"),
                       direction_threshold_edit_);

  QHBoxLayout *layout = new QHBoxLayout;
  layout->addLayout(layout_left);
  layout->addLayout(layout_right);

  gp->setLayout(layout);
  return gp;
}


} // namespace soax
