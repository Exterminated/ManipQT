#include "data_saver.h"

Data_Saver::Data_Saver()
{
    this->filename="";
    this->path="";
}

void Data_Saver::save_file(bool type){
    //false - txt
    //true - xml
    //TODO: придумать как реализовать сохранение файла, чтобы окно выбора всплывало только один раз
    if(type){

    }
    else{

    }
}
void Data_Saver::set_filename(QString filename){
    this->filename=filename;
}
QString Data_Saver::get_filename(){
    return filename;
}

