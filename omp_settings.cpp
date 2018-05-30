#include "omp_settings.h"
#include "ui_omp_settings.h"

omp_settings::omp_settings(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::omp_settings)
{
    ui->setupUi(this);
}

omp_settings::~omp_settings()
{
    delete ui;
}
