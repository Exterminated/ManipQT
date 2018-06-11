// заголовочный файл класса ManipCalculations.h
// интерфейс класса

#ifndef MANIPCALCULATIONS_H
#define MANIPCALCULATIONS_H

#define _USE_MATH_DEFINES


#include "cmath"
#include "omp_settings.h"
#include <conio.h>
#include <iostream>
#include <QDebug>
#include <omp.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include "alglib\optimization.h"

using namespace std;
using namespace alglib;

class ManipCalculations
{
  private:
  /* список свойств и методов для использования внутри класса */
  void function1_fvec(const real_1d_array &x, real_1d_array &fi, void *ptr);

  double xb, y_0, za, xc, yb, zd, O1A, OA1, l1, l2, l3, OO1;
  double  fi_angle, l4, A, gamma1, gamma0, gamma, xm0, ym0, zm0, xe0, ye0, ze0, psi1, psi2, alpha;
  double xmk, ymk, zmk, xek, yek, zek;
  int alpha_0, alpha_23, a, alpha_13, alpha_33, b;
  double OK, OA, OB, DK, Imin;

  double RungeKutta4(double x, double y,double h,double f(double,double));
  double f(double fi);

  double lk[4];
  double Time[4];
  double XP[4];
  double YP[4];
  double ZP[4];

  double Vmax = 25.0;
  int ptrsettings[3];

  public:
  /* список методов доступных другим функциям и объектам программы */
  ManipCalculations();
  ManipCalculations(double xb, double y_0,double za,double xc,double yb,
                    double zd, double O1A, double OA1, double l1, double l2, double l3,
                    double OO1,double fi_angle,double l4, double alpha, int alpha_0, int alpha_23,
                    int alpha_13,int alpha_33, int b);
  void calculatuons();
  void openmp_calculations();
  void setParams(double xb, double y_0,double za,double xc,double yb,
                 double zd, double O1A, double OA1, double l1, double l2, double l3,
                 double OO1,double fi_angle,double l4, double alpha, int alpha_0, int alpha_23,
                 int alpha_13,int alpha_33, int b);
  void set_lk(double l1k, double l2k, double l3k, double l4k);
  double *get_lk();

  void quickSortR(double* a, int N);
  double findMax(double* a, int N);

  void set_omp_settings(int num,int dynamic, int deapth);

  double *ptrTime = Time;
  double *ptrXP = XP;
  double *ptrYP = YP;
  double *ptrZP = ZP;

  double calculatuons_time = 0.0;

  protected:
  /*список средств, доступных при наследовании*/
};

#endif // MANIPCALCULATIONS_H
