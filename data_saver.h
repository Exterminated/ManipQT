#ifndef DATA_SAVER_H
#define DATA_SAVER_H

#include <QString>

#include<QFileDialog>
#include <QXmlStreamWriter>
#include <QXmlStreamReader>
#include <QXmlStreamAttribute>
#include <QMessageBox>
#include <QFile>

#include "manipcalculations.h"

class Data_Saver
{
public:
    Data_Saver();
    QString get_data(bool type, ManipCalculations item, QFile file);
    void set_filename(QString filename);
    QString get_filename();
private:
    QString path;
    QString filename;
};

#endif // DATA_SAVER_H
