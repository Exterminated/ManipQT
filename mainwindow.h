#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include <QDialog>
#include <QWidget>
#include <QMessageBox>
#include <QTime>

#include "data_saver.h"
#include "manipcalculations.h"
#include "omp_settings.h"
#include "about.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    omp_settings settings;
    about about_window;

private slots:
    void on_Calculate_pushButton_clicked();

    void on_actionTXT_triggered();

    void on_actionXML_triggered();

    void on_action_OpenMP_triggered();

    void on_action_about_triggered();

private:
    Ui::MainWindow *ui;
    void menue_visible(bool flag);
    void file_save(bool type);
    ManipCalculations manipcalculations;

    Data_Saver data_saver;

};
#endif // MAINWINDOW_H
