#ifndef DATA_SAVER_H
#define DATA_SAVER_H

#include <QString>

class Data_Saver
{
public:
    Data_Saver();
    void save_file(bool type);
    void set_filename(QString filename);
    QString get_filename();
private:
    QString path;
    QString filename;
};

#endif // DATA_SAVER_H
