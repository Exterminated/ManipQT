#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "manipcalculations.h"
#include <QMessageBox>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
//    QMessageBox msgBox;
//    msgBox.setWindowTitle("MessageBox Title");
//    msgBox.setText("You Clicked "+ ((QPushButton*)sender())->text());
//    msgBox.exec();
    ManipCalculations manipcalculations(
                //TODO порядок инициализации

                ui->alpha_spinBox,
                ui->alpha0_spinBox,
                ui->alpha13_spinBox,
                ui->alpha23_spinBox,
                ui->alpha33_spinBox,
                ui->b_spinBox,
                ui-> Xb_spinBox,
                ui-> Y0_spinBox,
                ui-> Za_spinBox,
                ui-> Xc_spinBox,
                ui-> Yb_spinBox,
                ui-> Zd_spinBox,
                ui-> O1A_spinBox,
                ui-> OA1_spinBox,
                ui-> l1_spinBox,
                ui-> l2_spinBox,
                ui-> l3_spinBox,
                ui-> O01_spinBox,
                ui-> fi_doubleSpinBox,
                ui-> l4_spinBox);
    manipcalculations.calculatuons();
}
