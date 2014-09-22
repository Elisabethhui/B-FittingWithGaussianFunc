#include"./Math.hpp"
#include<cstdio>

void Euler_Maclaurin(int Rm, int Nr, double * rm, double * w_rm){
//get Euler_maclaurin grid
	for(int n=1;n<=Nr;++n){
		rm[n-1]=(Rm*n*n)/((Nr+1.0-n)*(Nr+1.0-n));//parenthesis is necessry to convert expression to dboule precision
		w_rm[n-1]=2*(Rm*Rm*Rm)*(Nr+1.0)*n*n*n*n*n/std::pow(Nr+1.0-n, 7);
	}

}

std::complex<double> integral_3d_infinit(std::complex<double> (*eta)(double x,double y,double z,void * par), std::complex<double> (*f)(double x,double y,double z, void * par),void * par1, void * par2, int order, int Nr, double Rm){
//calculate the integration <eta|f>
//via lebedev integral and quadrature
//order is the lowest lebedev order
//Nr is the number of point for Euler-Maclaurin radius
for(int n=0;n<65;++n){
	if(available_table(n)>0 && order_table(n)>order ){
			order=order_table(n);
			break;
	}
}

double x[order], y[order],z[order],w[order];
ld_by_order(order,x,y,z,w);
double r[Nr],w_r[Nr];
Euler_Maclaurin(Rm,Nr,r,w_r);
std::complex<double> rslt(0);
std::complex<double> eta_val, f_val;
double xyz[3];
for(int n=0;n<Nr;++n){
	for(int m=0;m<order;++m){
		xyz[0]=x[m]*r[n];
		xyz[1]=y[m]*r[n];
		xyz[2]=z[m]*r[n];
		eta_val=eta(xyz[0],xyz[1],xyz[2], par1);	
		f_val=f(xyz[0],xyz[1],xyz[2], par2);
		rslt+=std::conj(eta_val)*f_val*w[m]*w_r[n];
//		printf("eta_val,w[m],w_r[n],rslt=%f,%f,%f,%f, xyz=%f,%f,%f\n",eta_val.real(),w[m],w_r[n],rslt.real(),xyz[0],xyz[1],xyz[2]);
	}
}

rslt*=4.0*3.14159265359;
return rslt;
}
