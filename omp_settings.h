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
    int num_threads;
    int deapth;
    int dynamic;

private slots:
    void on_buttonBox_accepted();

    void on_buttonBox_rejected();

private:
    Ui::omp_settings *ui;
    void set_saved_settings();
};

#endif // OMP_SETTINGS_H
