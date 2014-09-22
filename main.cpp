#include "./Math/Math.hpp"
#include<iostream>
#include<cstdio>
#include<cmath>
#include"IO/IO.hpp"
#include"Parallelization/Parallelization.hpp"
#include"./GaussFitting/GaussFitting.hpp"

std::complex<double>  fn(double x, double y, double z, void * par){
	return 1.0;
}

int main(int argc, char * argv[]){
	IO io;
	Parallelization parallel;
	parallel.Init(argc, argv);
	std::string file_name="inp.txt";
	io.Init(file_name, &parallel);
	GaussFitting fit;
	void * fn_par(NULL);
	fit.Init(& io, & fn, & fn_par);
	fit.Calc_g_g_mat();
	fit.Calc_g_f_mat();
	fit.OutputMat();
	return 0;
}
