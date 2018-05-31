#include "omp_settings.h"
#include "ui_omp_settings.h"

omp_settings::omp_settings(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::omp_settings)
{
    ui->setupUi(this);
//    set_saved_settings();
}

//omp_settings::show(){
//    set_saved_settings();
//    this->show();
//}

omp_settings::~omp_settings()
{
    delete ui;
}

void omp_settings::set_saved_settings(){
    if(num_threads!=ui->num_threads_spinBox->value()) ui->num_threads_spinBox->setValue(num_threads);
    if(deapth!=ui->lvl_spinBox->value())ui->lvl_spinBox->setValue(deapth);
    if(dynamic>0)ui->dynamic_checkBox->setChecked(true);
    else ui->dynamic_checkBox->setChecked(false);
}

void omp_settings::on_buttonBox_accepted()
{
    num_threads=ui->num_threads_spinBox->value();
    deapth=ui->lvl_spinBox->value();
    if(ui->dynamic_checkBox->isChecked()) dynamic=1;
    else dynamic=0;
}

void omp_settings::on_buttonBox_rejected()
{
    this->close();
}
