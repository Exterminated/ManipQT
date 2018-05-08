// файл реализации класса ManipCalculations.cpp
#include "manipcalculations.h"

// подключаем интерфейс класса к файлу его реализации
ManipCalculations::ManipCalculations(){
    qDebug()<<"New empty ManipCalculations object";
}
ManipCalculations::ManipCalculations(double xb, double y_0, double za, double xc, double yb,  double zd, double O1A, double OA1, double l1, double l2, double l3, double OO1, double fi_angle, double l4,  double alpha, int alpha_0, int alpha_23, int alpha_13, int alpha_33, int b){
    qDebug()<<"New ManipCalculations object";
    setParams(xb, y_0, za, xc, yb, zd, O1A, OA1, l1, l2, l3, OO1, fi_angle, l4, alpha, alpha_0, alpha_23, alpha_13, alpha_33, b);

}

void ManipCalculations::setParams(double xb, double y_0, double za, double xc, double yb, double zd, double O1A, double OA1, double l1, double l2, double l3, double OO1, double fi_angle, double l4, double alpha, int alpha_0, int alpha_23, int alpha_13, int alpha_33, int b){

    //this->a=a;
    this->alpha=alpha;
    this->alpha_0=alpha_0;
    this->alpha_13=alpha_13;
    this->alpha_23=alpha_23;
    this->alpha_33=alpha_33;
    this->b=b;
    this-> xb = xb;
    this-> y_0 =y_0;
    this-> za = za;
    this-> xc = xc;
    this-> yb=yb;
    this-> zd=zd;
    this-> O1A=O1A;
    this-> OA1=OA1;
    this-> l1=l1;
    this-> l2=l2;
    this-> l3=l3;
    this-> OO1=OO1;
    this-> fi_angle=fi_angle;
    this-> l4=l4;

    }

void ManipCalculations::function1_fvec(const real_1d_array &x, real_1d_array &fi, void *ptr)
{
    // this callback calculates

    fi[0] = pow((sqrt((xmk*xmk + (ymk + OA1*sin(fi_angle))*(ymk + OA1*sin(fi_angle)) + ((zmk - OA1*cos(fi_angle))*(zmk - OA1*cos(fi_angle))))) - l1), 2.0);
    fi[1] = pow((sqrt((OK - OA1*sin(fi_angle))*(OK - OA1*sin(fi_angle)) + (OA1*cos(fi_angle) - DK)*(OA1*cos(fi_angle) - DK)) - l4), 2.0);

    //fi[0] = 10 * pow(x[0] + 3, 2);
    //fi[1] = pow(x[1] - 3, 2);

}

void ManipCalculations::openmp_calculations(){}

void ManipCalculations::calculatuons() {
    real_1d_array x = "[0,0]";
    //double epsg = 0.0000000001;
    double epsg = 0.002;
    double epsf = 0;
    double epsx = 0;
    ae_int_t maxits = 0;
    minlmstate state;
    minlmreport rep;

    //Исходные конструктивные данные
//    xb = 355.0;
//    y_0 = 755.0;
//    za = 750.0;
//    xc = -355.0;
//    yb = y_0;
//    zd = -40.0;
//    O1A = za;
//    OA1 = 750.0;
//    l1 = 1400.0;
//    l2 = 1500.0;
//    l3 = 1352.0;
//    fi_angle = 0.323;
//    OO1 = zd;

    //Начальные координаты точки М и длина 4 звена
    l4 = sqrt((za*sin(M_PI_2 - fi_angle) - OO1)*(za*sin(M_PI_2 - fi_angle) - OO1) + (yb - za*cos(M_PI_2 - fi_angle))*(yb - za*cos(M_PI_2 - fi_angle)));
    A = 0.5*l2*l2 + 0.5*l3*l3 - l1*l1;
    xm0 = (l3*l3 - l2*l2) / (4 * xb);
    ym0 = sqrt((l1*l1) - (((l3*l3 - l2*l2)*(l3*l3 - l2*l2)) / (16 * xb * xb)) - (((A - xb*xb - za*za)*(A - xb*xb - za*za)) / (4 * za*za)))*cos(fi_angle) - ((A - xb*xb + za*za) / (2 * za))*sin(fi_angle);
    zm0 = sqrt((l1*l1) - (((l3*l3 - l2*l2)*(l3*l3 - l2*l2)) / (16 * xb * xb)) - (((A - xb*xb - za*za)*(A - xb*xb - za*za)) / (4 * za*za)))*sin(fi_angle) + ((A - xb*xb + za*za) / (2 * za))*cos(fi_angle);
    gamma0 = atan((-1 * xm0) / (ym0 + OA1*sin(fi_angle)));
    qDebug() << "l4: " << l4 ;
    qDebug() << "A: " << A ;
    qDebug() << "xm0: " << xm0 ;
    qDebug() << "ym0: " << ym0 ;
    qDebug() << "zm0: " << zm0 ;
    qDebug() << "gamma0: " << gamma0 ;
    //Начальные координаты захвата E, направляющих косинусов, зависящих от тех.процесса
//    alpha_0 = 0;
//    alpha_23 = 1;
//    a = 90;
//    alpha_13 = 0;
//    alpha_33 = 0;
//    b = 220;

    xe0 = xm0 + b*alpha_13 - a*cos(alpha_0)*sin(gamma0);
    ye0 = ym0 + b*alpha_23 + a*cos(alpha_0)*cos(gamma0);
    ze0 = zm0 + b*alpha_33 + a*sin(alpha_0);

    qDebug() << "xe0: " << xe0 ;
    qDebug() << "ye0: " << ye0 ;
    qDebug() << "ze0: " << ze0 ;

    //Конечные координаты захвата Е
    gamma = gamma0;
    xek = 300;
    yek = 1500;
    zek = -100;

    //Определяем углы psi1, gamma
    psi1 = asin(-alpha_23*sin(gamma) + alpha_13*cos(gamma));
    psi2 = asin(-alpha_23*sin(gamma) - alpha_13*cos(gamma));
    alpha = asin(alpha_33 / cos(psi1));
    gamma1 = atan((-xek) / yek);

    qDebug() << "psi1: " << psi1 ;
    qDebug() << "psi2: " << psi2 ;
    qDebug() << "alpha: " << alpha ;
    qDebug() << "gamma1: " << gamma1 ;

    xmk = xek - b*alpha_13 + a*cos(alpha)*sin(gamma1);
    ymk = yek - b*alpha_23 - a*cos(alpha)*cos(gamma1);
    zmk = zek - b*alpha_33 - a*sin(alpha);

    qDebug() << "xmk: " << xmk ;
    qDebug() << "ymk: " << ymk ;
    qDebug() << "zmk: " << zmk ;

    //Оптимизация угла фи

    OK = y_0;
    OA = za;
    OB = xb;
    DK = zd;
    OA1 = za;

//    minlmcreatev(2, x, 0.002, state);
//    minlmsetcond(state, epsg, epsf, epsx, maxits);
//    minlmoptimize(state, function1_fvec);
//    minlmresults(state, x, rep);

    //Imin = pow((sqrt((xmk*xmk + (ymk + OA1*sin(fi))*(ymk + OA1*sin(fi)) + ((zmk - OA1*cos(fi))*(zmk - OA1*cos(fi))))) - l1), 2.0) + pow((sqrt((OK - OA1*sin(fi))*(OK - OA1*sin(fi)) + (OA1*cos(fi) - DK)*(OA1*cos(fi) - DK)) - l4), 2.0);

    qDebug() <<"terminationtype: "<< int(rep.terminationtype);
    qDebug() <<"x: "<< x.tostring(2).c_str();

    _getch();
}
