#ifndef MATH_HPP_
#define MATH_HPP_
#include "./sphere_lebedev_rule.hpp"
#include<complex>

void Euler_Maclaurin(int Rm, int Nr, double * rm, double * w_rm);

std::complex<double> integral_3d_infinit(std::complex<double> (*eta)(double x,double y,double z, void * par), std::complex<double> (*f)(double x,double y,double z, void * par),void * par1, void * par2, int order, int Nr, double Rm);

#endif
