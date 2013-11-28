#include <QtGui>
#include "view_options_dialog.h"

namespace soax {

ViewOptionsDialog::ViewOptionsDialog(QWidget *parent) : QDialog(parent) {
  QVBoxLayout *layout = new QVBoxLayout;
  layout->addWidget(this->CreateSlicePlanesGroup());
  layout->addWidget(this->CreateMIPGroup());
  button_box_ = new QDialogButtonBox(QDialogButtonBox::Ok |
                                     QDialogButtonBox::Cancel);

  layout->addWidget(button_box_);
  setLayout(layout);
  setWindowTitle(tr("View Options"));
  connect(button_box_, SIGNAL(accepted()), this, SLOT(accept()));
  connect(button_box_, SIGNAL(rejected()), this, SLOT(reject()));
}

QGroupBox * ViewOptionsDialog::CreateSlicePlanesGroup() {
  QGroupBox *gb = new QGroupBox(tr("Window/Level for Slice Planes"));
  window_edit_ = new QLineEdit("0");
  level_edit_ = new QLineEdit("0");

  connect(window_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(level_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));

  QLabel *window_label = new QLabel(tr("Win"));
  QLabel *level_label = new QLabel(tr("Lev"));

  QHBoxLayout *layout = new QHBoxLayout;
  layout->addWidget(window_label);
  layout->addWidget(window_edit_);
  layout->addWidget(level_label);
  layout->addWidget(level_edit_);

  gb->setLayout(layout);
  return gb;
}

QGroupBox * ViewOptionsDialog::CreateMIPGroup() {
  QGroupBox *gb = new QGroupBox(tr("Intensity Range for MIP Rendering"));
  min_intensity_edit_ = new QLineEdit("0");
  max_intensity_edit_ = new QLineEdit("0");

  connect(min_intensity_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));
  connect(max_intensity_edit_, SIGNAL(textChanged(const QString &)),
          this, SLOT(EnableOKButton()));

  QLabel *min_intensity_label = new QLabel(tr("Min"));
  QLabel *max_intensity_label = new QLabel(tr("Max"));

  QHBoxLayout *layout = new QHBoxLayout;
  layout->addWidget(min_intensity_label);
  layout->addWidget(min_intensity_edit_);
  layout->addWidget(max_intensity_label);
  layout->addWidget(max_intensity_edit_);
  //layout->addStretch();
  gb->setLayout(layout);
  return gb;
}

void ViewOptionsDialog::EnableOKButton() {
  button_box_->button(QDialogButtonBox::Ok)->setEnabled(true);
}

void ViewOptionsDialog::DisableOKButton() {
  button_box_->button(QDialogButtonBox::Ok)->setEnabled(false);
}


double ViewOptionsDialog::GetWindow() {
  return window_edit_->text().toDouble();
}

double ViewOptionsDialog::GetLevel() {
  return level_edit_->text().toDouble();
}

void ViewOptionsDialog::SetWindow(double window) {
  QString s;
  s.setNum(window);
  window_edit_->setText(s);
}

void ViewOptionsDialog::SetLevel(double level) {
  QString s;
  s.setNum(level);
  level_edit_->setText(s);
}

double ViewOptionsDialog::GetMinIntensity() {
  return min_intensity_edit_->text().toDouble();
}

double ViewOptionsDialog::GetMaxIntensity() {
  return max_intensity_edit_->text().toDouble();
}

void ViewOptionsDialog::SetMinIntensity(double min) {
  QString s;
  s.setNum(min);
  min_intensity_edit_->setText(s);
}

void ViewOptionsDialog::SetMaxIntensity(double max) {
  QString s;
  s.setNum(max);
  max_intensity_edit_->setText(s);
}

} // namespace soax
