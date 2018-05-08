#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include<QFileDialog>

#include "data_saver.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_Calculate_pushButton_clicked();

    void on_actionTXT_triggered();

    void on_actionXML_triggered();

private:
    Ui::MainWindow *ui;
    void menue_visible(bool flag);
    void file_save();

    Data_Saver data_saver;

};

#endif // MAINWINDOW_H
