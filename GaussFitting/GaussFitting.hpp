#ifndef GAUSSFITTING_HPP
#define GAUSSFITTING_HPP
#include "../Math/Math.hpp"
#include<iostream>
#include<cstdio>
#include<cmath>
#include<complex>
#include"../IO/IO.hpp"

class GaussFitting{
public:
	typedef struct{ //Primitive Gaussian Basis.
		std::complex<double> coeff;
		std::complex<double> alpha;
		std::complex<double> R[3];
		int					 L[3];
		std::complex<double> norm;
	} PrimitiveGauss;
	typedef std::vector<PrimitiveGauss> ContractGauss;

	std::vector<ContractGauss>			basis_;
	std::vector<std::complex<double> >  g_f_mat_; //matrix <g|f>
	std::vector<std::complex<double> >  g_g_mat_; //matrix <g|g>
	IO									* io;
	std::complex<double> 				(*fn_)(double x, double y, double z, void * par);
	void								* fn_par_;
public:  
	void			ParseGaussianBasis();
	void 			FittingMatrixGenerate();
	int				Init(IO * io, std::complex<double>(*pfn)(double, double, double, void *), void * fn_par);
	void 			AddNewBasis(const std::vector<int> & L,const std::vector<std::complex<double> > & alpha, const std::vector<std::vector<double> > & coeff, double R[3]);
	int	 			InputParameter();
	void			CalcPrimitiveGaussNorm();
	void			AdjustBasisSize(const std::vector<int> & L);
	void			Calc_g_f_mat();
	void			Calc_g_g_mat();
	void			OutputMat();
};

#endif
