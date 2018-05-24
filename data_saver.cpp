#include "data_saver.h"

Data_Saver::Data_Saver()
{
    this->filename="";
    this->path="";
}

//QString Data_Saver::save_file(bool type, ManipCalculations item, QFile file){
//    //false - txt
//    //true - xml

//    //TODO: придумать как реализовать сохранение файла, чтобы окно выбора всплывало только один раз
////    if(type){
////    }
////    else{

////        for(int i=0;i<4;i++){
////            outString+=QString::number(item.Time[i])+";";
////        }
////    }
////    return outString;
//}
void Data_Saver::set_filename(QString filename){
    this->filename=filename;
}
QString Data_Saver::get_filename(){
    return filename;
}

