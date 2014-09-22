#include "../Math/Math.hpp"
#include<iostream>
#include<cstdio>
#include<cmath>
#include<complex>
#include"../IO/IO.hpp"
#include"./GaussFitting.hpp"
#include<iomanip>

std::complex<double> GaussianFn(double x, double y, double z, void * par){
	GaussFitting::PrimitiveGauss * parameter=static_cast<GaussFitting::PrimitiveGauss *>(par);
	std::complex<double> rR[3];
	rR[0]=x-parameter->R[0];
	rR[1]=y-parameter->R[1];
	rR[2]=z-parameter->R[2];
	return std::pow(rR[0],parameter->L[0])*std::pow(rR[1],parameter->L[1])*std::pow(rR[2],parameter->L[2])*\
			exp(-parameter->alpha*(rR[0]*rR[0]+rR[1]*rR[1]+rR[2]*rR[2]));
}

int GaussFitting::InputParameter(){
	ParseGaussianBasis();
	return 0;
}

void GaussFitting::ParseGaussianBasis(){
//This function is to parse the basis block GaussianBasis
	std::string varName="GaussianBasis";
	if(!io->inpParser.IfVariableExist(varName)){
		io->PPrintf("Gaussian basis set has to be provided\n");		
		exit(1);
	}
	int row,col, ri;
	io->inpParser.GetBlockRowLen(varName, row);	
	std::string 	str;
	double			R[3];
	bool			grpFlag=false;//if reading gaussian basis with same center
	std::vector<ContractGauss>	& basis(basis_);
	ri=0;
	io->inpParser.GetBlockValue(varName, ri+1,1,"string",&str);
	if(str.compare("R")!=0){
		io->PPrintf("Gaussianbasis: row %d char %s is incorrect\n",ri+1, str.c_str());
		exit(1);
	}
	while(1){
		if(!grpFlag){//reading the center first
			grpFlag=true;
			io->inpParser.GetBlockColLen(varName, ri+1, col);
			if(col!=4){
				io->PPrintf("row %d format is incorrect\n", ri+1);
				exit(1);
			}
			for(int ci=0;ci<3;++ci){
				io->inpParser.GetBlockValue(varName, ri+1,ci+2,"double",&R[ci]);	
			}
			ri++;
		}else{//read contract Gaussian basis;
			io->inpParser.GetBlockValue(varName, ri+1,1,"string",&str);
			if(str.size()==0){
				io->PPrintf("Gaussian basis invalid format at row %d\n", ri+1);
				exit(1);
			}
			std::vector<int> L;
			if(str[0]=='L'){//if it is L, momentum is determined by the next row column
				if(str.size()>1){//L should not be combined with any other momentum
					io->PPrintf("Gaussian basis invalid format at row %d\n", ri+1);
					exit(1);
				}
				io->inpParser.GetBlockColLen(varName,ri+2,col);	
				L.resize(col-1);
				for(int n=0;n<L.size();++n)
					L[n]=n;
			}else{//if not L, it is should be the combination of S, P, D, F, G, H, other higher momentum is not supported;
				L.resize(str.size());
				for(int n=0;n<str.size();++n){
					switch(str[n]){
					case 'S':
						L[n]=0;
						break;
					case 'P':
						L[n]=1;
						break;
					case 'D':
						L[n]=2;
						break;
					case 'F':
						L[n]=3;
						break;
					case 'G':
						L[n]=4;
						break;
					case 'H':
						L[n]=5;
						break;
					default:
						io->PPrintf("in GaussianBasis invalid momentum format\n");
						exit(1);
					}
				}
			}

			int cn;//primitive number;
			io->inpParser.GetBlockValue(varName, ri+1,2,"int", &cn);	
			ri++;
			std::vector<std::complex<double> > alpha(cn);
			std::vector<std::vector<double> > coeff;
			coeff.resize(cn);
			for(int n=0;n<cn;++n){
				coeff[n].resize(L.size());
			}

			//now read the coefficient;
			for(int n=0;n<cn;++n){//cn is the coefficient number;
				io->inpParser.GetBlockValue(varName, ri+1,1,"complex",&alpha[n]);
				for(int m=0;m<L.size();++m){
					io->inpParser.GetBlockValue(varName, ri+1,m+2,"double",&coeff[n][m]);//get coefficient
				}
				ri++;
			}
			//now add new basis 
			AddNewBasis(L, alpha, coeff, R);
			//now check if move to a new group
			if(ri>=row){
				break;
			}else{
				io->inpParser.GetBlockValue(varName, ri+1,1,"string",&str);
				grpFlag=(str[0]=='R')? false:true;
			}
		}
	}
	if(basis.size()==0){
		io->PPrintf("Gaussian basis set cannot be empty\n");
		exit(1);
	}
}

void GaussFitting::AddNewBasis(const std::vector<int> & L,const std::vector<std::complex<double> > & alpha, const std::vector<std::vector<double> > & coeff, double R[3]){
	std::vector<ContractGauss>	&basis(basis_);
	int idx=basis.size();
	PrimitiveGauss *p;
	AdjustBasisSize(L);
	for(int n=0;n<L.size();++n){
		for(int lx=0;lx<=L[n];++lx){
			for(int ly=0;ly<=L[n];++ly){
				for(int lz=0;lz<=L[n];++lz){
					if(lx+ly+lz==L[n]){
						basis[idx].resize(alpha.size());
						for(int m=0;m<alpha.size();++m){
							p=&basis[idx][m];
							p->alpha=alpha[m];
							p->coeff=coeff[m][n];
							p->R[0]=R[0];
							p->R[1]=R[1];
							p->R[2]=R[2];
							p->L[0]=lx;
							p->L[1]=ly;
							p->L[2]=lz;
						}
						++idx;
					}
				}
			}
		}
	}		
}

void GaussFitting::CalcPrimitiveGaussNorm(){
//after reading Gaussian basis set, calc the norm of each primitive gaussian basis
	int * L;	
	for(int n=0;n<basis_.size();++n){
		for(int m=0;m<basis_[n].size();++m){
			L=basis_[n][m].L;
			basis_[n][m].norm=(std::pow(2.0*basis_[n][m].alpha,(2.0*(L[0]+L[1]+L[2])+3.0)/4.0))/sqrt(tgamma(L[0]+0.5)*tgamma(L[1]+0.5)*tgamma(L[2]+0.5));			
		}
	}
}

void GaussFitting::AdjustBasisSize(const std::vector<int> & L){
	int num=0;
	for(int n=0;n<L.size();++n){
			num+=(L[n]==0)? 1:3*L[n];
	}
	basis_.resize(basis_.size()+num);
}

int GaussFitting::Init(IO * io, std::complex<double>(*pfn)(double, double, double, void *), void * fn_par){
	this->io=io;
	fn_=pfn;
	fn_par_=fn_par;
	if(InputParameter()) return 1;
	CalcPrimitiveGaussNorm();
	return 0;
}

void GaussFitting::Calc_g_f_mat(){
	int size=basis_.size();
	std::complex<double> rslt;
	g_f_mat_.resize(size,0);
	for(int n=0;n<basis_.size();++n){//for each contract basis
		rslt=0;
		for(int m=0;m<basis_[n].size();++m){//for each primitive basis
			rslt+=integral_3d_infinit(&GaussianFn,fn_,&basis_[n][m],fn_par_,100,100,1)*basis_[n][m].norm*basis_[n][m].coeff;//todo: order might need to change
		}
		g_f_mat_[n]=rslt;
	}	
}

void GaussFitting::Calc_g_g_mat(){
	int size=basis_.size();
	std::complex<double> rslt,temp;
	g_g_mat_.resize(size*size,0);
	for(int n=0;n<basis_.size();++n){//for each contract basis
		for(int k=n;k<basis_.size();++k){
			rslt=0;
			for(int m1=0;m1<basis_[n].size();++m1){//for each primitive basis
				for(int m2=0;m2<basis_[k].size();++m2){
					rslt+=integral_3d_infinit(&GaussianFn,&GaussianFn, &basis_[n][m1],&basis_[k][m2],100,100,10)*\
						  std::conj(basis_[n][m1].norm*basis_[n][m1].coeff)*basis_[k][m2].norm*basis_[k][m2].coeff;//todo: order might need to change
				}
			}
			g_g_mat_[n*size+k]=rslt;
			g_g_mat_[k*size+n]=std::conj(rslt);
		}
	}	
}

void GaussFitting::OutputMat(){
	std::string file_name="g_g_mat.txt";
	std::ofstream * ofile=io->outPut.AddOutputFile(file_name);
	int size=basis_.size();		
	int prec(20),width(30);
	(*ofile)<<"g_g_mat: "<<size<<"x"<<size<<std::endl;
	for(int n=0;n<size;++n){
		for(int c=0;c<size;++c){
			(*ofile)<<std::setprecision(prec)<<std::setw(width)<<g_g_mat_[n*size+c];
		}
		(*ofile)<<std::endl;
	}
	io->outPut.CloseFile(file_name);
	
	file_name="g_f_mat.txt";
	ofile=io->outPut.AddOutputFile(file_name);
	(*ofile)<<"g_f_mat: "<<size<<std::endl;
	for(int n=0;n<size;++n){
		(*ofile)<<std::setprecision(prec)<<std::setw(width)<<g_f_mat_[n]<<std::endl;
	}
	io->outPut.CloseFile(file_name);
}


