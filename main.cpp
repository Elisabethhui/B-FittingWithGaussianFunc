#include "./Math/Math.hpp"
#include<iostream>
#include<cstdio>
#include<cmath>
#include<complex>
struct GaussPar{
	int l[3];
	std::complex<double> R[3];
	std::complex<double> alpha;
	GaussPar();
};

GaussPar::GaussPar(){
	l[0]=0;
	l[1]=0;
	l[2]=0;
	R[0]=0;
	R[1]=0;
	R[2]=0;
	alpha=0;
}




std::complex<double> GaussianFn(double x, double y, double z, void * par){
	GaussPar * parameter=static_cast<GaussPar *>(par);
	std::complex<double> rR[3];
	rR[0]=x-parameter->R[0];
	rR[1]=x-parameter->R[1];
	rR[2]=x-parameter->R[2];
	return std::pow(rR[0],parameter->l[0])*std::pow(rR[1],parameter->l[1])*std::pow(rR[2],parameter->l[2])*\
			exp(-parameter->alpha*(rR[0]*rR[0]+rR[1]*rR[1]+rR[2]*rR[2]));
}

int main(){
	std::complex<double> rslt;
	GaussPar par1;
	par1.alpha=3.0;
	rslt=integral_3d_infinit(&GaussianFn,&GaussianFn,&par1, & par1, 100, 100,1);
	printf("the integration of Gaussian function is (%f,%f)\n",rslt.real(),rslt.imag());
	return 0;
}
