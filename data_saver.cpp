#include "data_saver.h"

Data_Saver::Data_Saver()
{
    this->filename="";
    this->path="";
}
void Data_Saver::set_filename(QString filename){
    this->filename=filename;
}
QString Data_Saver::get_filename(){
    return filename;
}
QString Data_Saver::get_full_results_txt(ManipCalculations item){
    QString results=QString::null;
    for(int i=0;i<4;i++)results+="X["+QString::number(i)+"]:"+QString::number(item.ptrXP[i])+"; ";
    results+="\n";
    for(int i=0;i<4;i++)results+="Y["+QString::number(i)+"]:"+QString::number(item.ptrYP[i])+"; ";
    results+="\n";
    for(int i=0;i<4;i++)results+="Z["+QString::number(i)+"]:"+QString::number(item.ptrZP[i])+"; ";
    results+="\n";
    for(int i=0;i<4;i++)results+="T["+QString::number(i)+"]:"+QString::number(item.ptrTime[i])+"; ";
    results+="\n";
    return results;
}

