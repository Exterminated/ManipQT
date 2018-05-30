#ifndef OMP_SETTINGS_H
#define OMP_SETTINGS_H

#include <QDialog>

namespace Ui {
class omp_settings;
}

class omp_settings : public QDialog
{
    Q_OBJECT

public:
    explicit omp_settings(QWidget *parent = 0);
    ~omp_settings();

private:
    Ui::omp_settings *ui;
};

#endif // OMP_SETTINGS_H
