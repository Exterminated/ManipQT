#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "manipcalculations.h"
#include <QMessageBox>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    menue_visible(false);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::menue_visible(bool flag){
     ui->Export_menu->setEnabled(flag);
     ui->actionTXT->setEnabled(flag);
     ui->actionXML->setEnabled(flag);
}

void MainWindow::on_Calculate_pushButton_clicked()
{
    ManipCalculations manipcalculations(
                ui-> Xb_spinBox->value(),
                ui-> Y0_spinBox->value(),
                ui-> Za_spinBox->value(),
                ui-> Xc_spinBox->value(),
                ui-> Yb_spinBox->value(),
                ui-> Zd_spinBox->value(),
                ui-> O1A_spinBox->value(),
                ui-> OA1_spinBox->value(),
                ui-> l1_spinBox->value(),
                ui-> l2_spinBox->value(),
                ui-> l3_spinBox->value(),
                ui-> O01_spinBox->value(),
                ui-> fi_doubleSpinBox->value(),
                ui-> l4_spinBox->value(),
                ui->alpha_spinBox->value(),
                ui->alpha0_spinBox->value(),
                ui->alpha23_spinBox->value(),
                ui->alpha13_spinBox->value(),
                ui->alpha33_spinBox->value(),
                ui->b_spinBox->value());
    manipcalculations.calculatuons();
    menue_visible(true);
}

void MainWindow::on_actionTXT_triggered()
{
    if(data_saver.get_filename().isEmpty()){
    data_saver.set_filename(QFileDialog::getSaveFileName(this,"File save","","Text file|txt")+".txt");
    }
    file_save();
}

void MainWindow::on_actionXML_triggered()
{
    if(data_saver.get_filename().isEmpty()){
    data_saver.set_filename(QFileDialog::getSaveFileName(this,"File save","","XML files|xml")+".xml");
    }
    file_save();
}

void MainWindow::file_save(){
    if(!data_saver.get_filename().isEmpty()){
        QFile file(data_saver.get_filename());
        if(!file.open(QIODevice::WriteOnly)){
            QMessageBox::information(this,"Unable to open file",file.errorString());
            return;
        }
        QDataStream out(&file);
        out.setVersion(QDataStream::Qt_5_10);
        //Получить данные из data_saver
        out<<"Test strin";
    }
    else{
        QMessageBox::information(this,"Unable to open file","File name dose not set");
    }

}
