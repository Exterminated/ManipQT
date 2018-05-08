// заголовочный файл класса ManipCalculations.h
// интерфейс класса

#ifndef MANIPCALCULATIONS_H
#define MANIPCALCULATIONS_H

#define _USE_MATH_DEFINES

#include "alglib\optimization.h"
#include "cmath"
#include <conio.h>
#include <iostream>
#include <QDebug>

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

  public:
  /* список методов доступных другим функциям и объектам программы */
  ManipCalculations();
  ManipCalculations(double xb, double y_0,double za,double xc,double yb,
                    double zd, double O1A, double OA1, double l1, double l2, double l3,
                    double OO1,double fi_angle,double l4, double alpha, int alpha_0, int alpha_23,
                    int alpha_13,int alpha_33, int b);
  void calculatuons();
  void setParams(double xb, double y_0,double za,double xc,double yb,
                 double zd, double O1A, double OA1, double l1, double l2, double l3,
                 double OO1,double fi_angle,double l4, double alpha, int alpha_0, int alpha_23,
                 int alpha_13,int alpha_33, int b);

  protected:
  /*список средств, доступных при наследовании*/
};

#endif // MANIPCALCULATIONS_H
