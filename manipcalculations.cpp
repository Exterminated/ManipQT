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

void ManipCalculations::set_lk(double l1k, double l2k, double l3k, double l4k){
    this->lk[0]=l1k;
    this->lk[1]=l2k;
    this->lk[2]=l3k;
    this->lk[3]=l4k;
}
void ManipCalculations::set_omp_settings(int num,int dynamic,int deapth){
    this->ptrsettings[0]=num;
    this->ptrsettings[1]=dynamic;
    this->ptrsettings[2]=deapth;
}
void ManipCalculations::openmp_calculations(){
    omp_set_dynamic(ptrsettings[1]);
    omp_set_num_threads(ptrsettings[0]);
    omp_set_nested(ptrsettings[2]);

    double Sx [4] = {};
    double Sy [4] = {};
    double Sz [4] = {};

    double tau [4] = {};

    #if defined(_OPENMP)
        qWarning("Compiled by an OpenMP-compliant implementation.\n");
        qWarning("The result of omp_get_num_threads %i\n", omp_get_num_threads());
    #endif

#pragma omp parallel{}
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

  //fi_angle = RungeKutta4(0.0,0.0,1,(fi[0],fi[1]));
        qDebug()<<"fi angle"<<fi_angle;
        //expected 0.072
        qDebug()<<"function f"<<ManipCalculations::f(0.0,0.0);
        fi_angle=0.072;


        //рассчет lk
        double l1k = sqrt(xmk*xmk+pow(ymk+OA1*sin(fi_angle),2.0)+pow(zmk-OA1*cos(fi_angle),2.0));
        double l2k = sqrt((xmk-xb)*(xmk-xb)+ymk*ymk+zmk*zmk);
        double l3k = sqrt((xmk+xb)*(xmk*xb)+ymk*ymk+zmk*zmk);
        double l4k = sqrt(pow(OK-OA*sin(fi_angle),2.0)+pow(OA*cos(fi_angle)-DK,2.0));
        this->set_lk(l1k,l2k,l3k,l4k);

        //Рассчет времени
        ptrTime[0] = abs(l1k-l1)/Vmax;
        ptrTime[1] = abs(l2k-l2)/Vmax;
        ptrTime[2] = abs(l3k-l3)/Vmax;
        ptrTime[3] = abs(l4k-l4)/Vmax;
        double Tlk = ManipCalculations::findMax(ptrTime,4);

        qDebug()<<"Max time is "<<Tlk;

        //рассчет координат
        double dx = xmk-xm0;
        double dy = ymk-ym0;
        double dz = zmk-zm0;

        double k1 = dy/dx;
        double k2 = dz/dx;
        double k3 = dz/dy;

        qDebug()<< "k1 "<<k1;
        qDebug()<< "k2 "<<k2;
        qDebug()<<"k3 "<<k3;

        double kx = sqrt(1+k1*k1+k2*k2);
        double ky = sqrt(1+(1/k1)*(1/k1)+k3*k3);
        double kz = sqrt(1+(1/k2)*(1/k2)+(1/k3)*(1/k3));

        qDebug()<<"D = "<<((1/(kx*kx))+(1/(ky*ky))+(1/(kz*kz)))<<"; Should be 1";

        double Skx = sqrt(1+k1*k1+k2*k2)*(xmk-xm0);
        double Sky = sqrt(1+((1/(k1*k1))+(k3*k3)))*(ymk-ym0);
        double Skz = sqrt(1+((1/(k1*k1))+(1/(k3*k3))))*(zmk-zm0);

        qDebug()<<"Skx: "<<Skx;
        qDebug()<<"Sky: "<<Sky;
        qDebug()<<"Skz: "<<Skz;

        //Se: 0 - x; 1 - y; 2 - z;

    #pragma omp for{
        for(int i = 0; i<4;i++){
            tau[i]=(i+1)/Tlk;
            //qDebug<"tau = "<<tau[i];
        }
}
#pragma omp for{

        for(int i=0;i<4;i++){
            Sx[i]=Skx*(10.0*pow(tau[i],3.0)-15.0*pow(tau[i],4)+6.0*pow(tau[i],5.0));
            Sy[i]=Sky*(10.0*pow(tau[i],3.0)-15.0*pow(tau[i],4)+6.0*pow(tau[i],5.0));
            Sz[i]=Skz*(10.0*pow(tau[i],3.0)-15.0*pow(tau[i],4)+6.0*pow(tau[i],5.0));
            //---------------------------//
            ptrXP[i]=(Sx[i]+xm0*kx)/kx;
            ptrYP[i]=(Sy[i]+ym0*ky)/ky;
            ptrZP[i]=(Sz[i]+zm0*kz)/kz;

            qDebug()<<"X["<<i<<"]: "<<ptrXP[i];
            qDebug()<<"Y["<<i<<"]: "<<ptrYP[i];
            qDebug()<<"Z["<<i<<"]: "<<ptrZP[i];
        }
}

}

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

//    qDebug() <<"terminationtype: "<< int(rep.terminationtype);
//    qDebug() <<"x: "<< x.tostring(2).c_str();



    //fi_angle = RungeKutta4(0.0,0.0,1,(fi[0],fi[1]));
    qDebug()<<"fi angle"<<fi_angle;
    //expected 0.072
    qDebug()<<"function f"<<ManipCalculations::f(0.0,0.0);
    fi_angle=0.072;


    //рассчет lk
    double l1k = sqrt(xmk*xmk+pow(ymk+OA1*sin(fi_angle),2.0)+pow(zmk-OA1*cos(fi_angle),2.0));
    double l2k = sqrt((xmk-xb)*(xmk-xb)+ymk*ymk+zmk*zmk);
    double l3k = sqrt((xmk+xb)*(xmk*xb)+ymk*ymk+zmk*zmk);
    double l4k = sqrt(pow(OK-OA*sin(fi_angle),2.0)+pow(OA*cos(fi_angle)-DK,2.0));
    this->set_lk(l1k,l2k,l3k,l4k);

    //Рассчет времени   
    ptrTime[0] = abs(l1k-l1)/Vmax;
    ptrTime[1] = abs(l2k-l2)/Vmax;
    ptrTime[2] = abs(l3k-l3)/Vmax;
    ptrTime[3] = abs(l4k-l4)/Vmax;
    double Tlk = ManipCalculations::findMax(ptrTime,4);

    qDebug()<<"Max time is "<<Tlk;

    //рассчет координат
//    OA=za;
//    OB=xb;
//    OK=y_0;
//    DK=zd;
//    OA1=za;

    double dx = xmk-xm0;
    double dy = ymk-ym0;
    double dz = zmk-zm0;

    double k1 = dy/dx;
    double k2 = dz/dx;
    double k3 = dz/dy;

    qDebug()<< "k1 "<<k1;
    qDebug()<< "k2 "<<k2;
    qDebug()<<"k3 "<<k3;

    double kx = sqrt(1+k1*k1+k2*k2);
    double ky = sqrt(1+(1/k1)*(1/k1)+k3*k3);
    double kz = sqrt(1+(1/k2)*(1/k2)+(1/k3)*(1/k3));

    qDebug()<<"D = "<<((1/(kx*kx))+(1/(ky*ky))+(1/(kz*kz)))<<"; Should be 1";

    double Skx = sqrt(1+k1*k1+k2*k2)*(xmk-xm0);
    double Sky = sqrt(1+((1/(k1*k1))+(k3*k3)))*(ymk-ym0);
    double Skz = sqrt(1+((1/(k1*k1))+(1/(k3*k3))))*(zmk-zm0);

    qDebug()<<"Skx: "<<Skx;
    qDebug()<<"Sky: "<<Sky;
    qDebug()<<"Skz: "<<Skz;

    //Se: 0 - x; 1 - y; 2 - z;
    double Sx [4] = {};
    double Sy [4] = {};
    double Sz [4] = {};

//    double *xp = XP;
//    double *yp = YP;
//    double *zp = ZP;
    double tau [4] = {};

    for(int i = 0; i<4;i++){
        tau[i]=(i+1)/Tlk;
        //qDebug<"tau = "<<tau[i];
    }

    for(int i=0;i<4;i++){
        Sx[i]=Skx*(10.0*pow(tau[i],3.0)-15.0*pow(tau[i],4)+6.0*pow(tau[i],5.0));
        Sy[i]=Sky*(10.0*pow(tau[i],3.0)-15.0*pow(tau[i],4)+6.0*pow(tau[i],5.0));
        Sz[i]=Skz*(10.0*pow(tau[i],3.0)-15.0*pow(tau[i],4)+6.0*pow(tau[i],5.0));
        //---------------------------//
        ptrXP[i]=(Sx[i]+xm0*kx)/kx;
        ptrYP[i]=(Sy[i]+ym0*ky)/ky;
        ptrZP[i]=(Sz[i]+zm0*kz)/kz;

        qDebug()<<"X["<<i<<"]: "<<ptrXP[i];
        qDebug()<<"Y["<<i<<"]: "<<ptrYP[i];
        qDebug()<<"Z["<<i<<"]: "<<ptrZP[i];
    }
}
void ManipCalculations::quickSortR(double *a, int N){
    int i = 0, j= N-1;
    double temp, p;
    p=a[N>>1];
    do
    {
        while(a[i]<p) i++;
           while(a[i]>p)j--;

           if(i<=j)
        {
           temp=a[i];
           a[i]=a[j];
           a[j]=temp;
           i++;j--;
        }
    }while (i<=j);
    if(j>0) ManipCalculations::quickSortR(a,j);
    if(N>i)ManipCalculations::quickSortR(a+i,N-i);
}
double ManipCalculations::findMax(double* a, int N){
    double max = 0.0;
    for(int i=0;i<0;i++){
        if(a[i]>max) max=a[i];
    }
    return max;
}
void ManipCalculations::function1_fvec(const real_1d_array &x, real_1d_array &fi, void *ptr)
{
    // this callback calculates

    fi[0] = pow((sqrt((xmk*xmk + (ymk + OA1*sin(fi_angle))*(ymk + OA1*sin(fi_angle)) + ((zmk - OA1*cos(fi_angle))*(zmk - OA1*cos(fi_angle))))) - l1), 2.0);
    fi[1] = pow((sqrt((OK - OA1*sin(fi_angle))*(OK - OA1*sin(fi_angle)) + (OA1*cos(fi_angle) - DK)*(OA1*cos(fi_angle) - DK)) - l4), 2.0);

    //fi[0] = 10 * pow(x[0] + 3, 2);
    //fi[1] = pow(x[1] - 3, 2);

}
//***************************************************************

            //** Метод Рунге-Кутта 4-го порядка точности.
            //** Параметры: x,y - типа float,
            //**            h   - шаг,
            //**            f   - функция сорж. диф. уравнение.

double ManipCalculations::RungeKutta4(double x, double y,double h,double f(double,double))
{
  double k1,k2,k3,k4;

  k1=h*f(x,y);
  k2=h*f(x+h/2.0,y+k1/2.0);
  k3=h*f(x+h/2.0,y+k2/2.0);
  k4=h*f(x+h,y+k3);
  return(y+(k1+2*k2+2*k3+k4)/6);
}

double ManipCalculations::f(double x, double y){
//    double fi[2] = {
//        (double)pow((sqrt((xmk*xmk + (ymk + OA1*sin(fi_angle))*(ymk + OA1*sin(fi_angle)) + ((zmk - OA1*cos(fi_angle))*(zmk - OA1*cos(fi_angle))))) - l1), 2.0),
//        (double)pow((sqrt((OK - OA1*sin(fi_angle))*(OK - OA1*sin(fi_angle)) + (OA1*cos(fi_angle) - DK)*(OA1*cos(fi_angle) - DK)) - l4), 2.0)
//    };
    double denominator_1, numerator_1, denominator_2, numerator_2;
    denominator_1 = (2.0*OA1*cos(fi_angle)*(ymk+OA1*sin(fi_angle))+2.0*OA1*sin(fi_angle)*(zmk-OA1*cos(fi_angle)))*(l1-sqrt(pow(ymk+OA1*sin(fi_angle),2.0)+pow(zmk-OA1*cos(fi_angle),2.0)));
    numerator_1 = sqrt(pow(ymk+OA1*sin(fi_angle),2.0)+pow(zmk-OA1*cos(fi_angle),2.0)+xmk*xmk);
    denominator_2 = (l4-sqrt(pow(DK-OA1*cos(fi_angle),2.0)+pow(OK-OA1*sin(fi_angle),2.0)))*(2*OA1*sin(fi_angle)*(DK-OA1*cos(fi_angle))-2.0*OA1*cos(fi_angle)*(OK-OA1*sin(fi_angle)));
    numerator_2 = sqrt(pow(DK-OA1*cos(fi_angle),2.0)+pow(OK-OA1*sin(fi_angle),2.0));
    return (denominator_1/numerator_1)+(denominator_2/numerator_2);
}
