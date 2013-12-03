#ifndef SOAX_VIEW_OPTIONS_DIALOG_H_
#define SOAX_VIEW_OPTIONS_DIALOG_H_

#include <QDialog>
#include "global.h"

class QLineEdit;
class QGroupBox;
class QDialogButtonBox;

namespace soax {

class ViewOptionsDialog : public QDialog {
  Q_OBJECT

 public:
  ViewOptionsDialog(QWidget * parent = NULL);

  double GetWindow() const;
  double GetLevel() const;

  void SetWindow(double window);
  void SetLevel(double level);

  double GetMinIntensity() const;
  double GetMaxIntensity() const;

  void SetMinIntensity(double min);
  void SetMaxIntensity(double max);

  double GetClipSpan() const;
  void SetClipSpan(double span);

  unsigned GetColorSegmentStep() const;
  void SetColorSegmentStep(unsigned step);

 public slots:
  void EnableOKButton();
  void DisableOKButton();

 private:
  QGroupBox *CreateSlicePlanesGroup();
  QGroupBox *CreateMIPGroup();
  QGroupBox *CreateClipGroup();
  QGroupBox *CreateColorOrientationGroup();

  QLineEdit *window_edit_;
  QLineEdit *level_edit_;
  QLineEdit *min_intensity_edit_;
  QLineEdit *max_intensity_edit_;
  QLineEdit *clip_span_edit_;
  QLineEdit *color_segment_step_edit_;

  QDialogButtonBox *button_box_;

  DISALLOW_COPY_AND_ASSIGN(ViewOptionsDialog);
};

} // namespace soax
#endif // SOAX_VIEW_OPTIONS_DIALOG_H_
