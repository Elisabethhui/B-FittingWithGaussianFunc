#ifndef _GAUSSIANBASIS_CPP_
#define _GAUSSIANBASIS_CPP_
#include "ContractGaussianBasis.hpp"
#include<iostream>
#include<cstring>
#include<iomanip>
#include<gsl/gsl_integration.h>

ContractGaussianBasis::ContractGaussianBasis():
shellPairCutoff(1e-12),signifShellPairRelDist(20),GMMaxM(30)
{
}


inline bool ContractGaussianBasis::CheckVar(const std::string & var)
{
	if(io->inpParser.IfVariableExist(var))
	{
		io->inpParser.ReserveVar(var,"ContractGaussianBasis");
		return true;
	}else{
		return false;
	}
}

void ContractGaussianBasis::InputParameters()
{//get parameters from input file
	//SignifShellPairRelDist
	
	ParseGaussianBasis();	

	std::string varName="ShellPairCutoff";
	if(CheckVar(varName))
	{
		io->inpParser.GetVariableValue(varName,"double",&signifShellPairRelDist);
	}else
	{
		signifShellPairRelDist=20; //[a.u.] default value;
	}

	//shellPairCutoff
	varName="shellPairCutoff";
	if(CheckVar(varName))
	{
		io->inpParser.GetVariableValue(varName,"double",&shellPairCutoff);
	}else
	{
		shellPairCutoff=1e-15; //[a.u.] default value;
	}
	//roundBracket
	varName="RoundBracket";
	if(CheckVar(varName))
	{
		io->inpParser.GetVariableValue(varName,"bool",&roundBracket);
	}else
	{
		roundBracket=true;
	}
	if(false==roundBracket){
		std::cerr<<"complex conjugated integral is not support yet"<<std::endl;
		exit(1);
	}
}

int ContractGaussianBasis::Init(IO * io,Species * species,size_t dim){
	SpaceBasis::Init(io, species, dim);
	InputParameters();
	CalcPrimitiveGaussNorm();


	matrixSize=basis.size();//for tddft, the matrix should be the same as size, but for two electron system, it is not true
	size=basis.size();
	BasisPartition(); //must be called after setting matrixSize;

//	find the possible maximum m for [0]^m
	int L[3]={0};
	for(int n=0;n<basis.size();++n){
		for(int m=0;m<basis[n].size();++m){
			for(int i=0;i<3;++i)
				if(L[i]<basis[n][m].L[i])  L[i]=basis[n][m].L[i];
		}
	}
	GMMaxM=(L[0]+L[1]+L[2])*4+2;//?? is it correct??
	
	io->PPrintf("ContractGaussianBasis is initiated\n");
	return 0;
}



void ContractGaussianBasis::help()
{

}

bool ContractGaussianBasis::Check(){
	bool rslt=false;
	if(NULL==io) return false;
	if(0==this->size) return false;
	if(GMMaxM==-1) rslt=false;
	return rslt;
}


void ContractGaussianBasis::GetLocalSpaceTermPiecen(int piecen, const std::string & type,std::vector<std::complex<double> > & rslt){
//has to rewrite; spaceBasis_abcd no longer usefull
	//if it is already calcualted, if not, calculate it.
	switch(elementType.at(type)){
	case 0: //passed test
		CalcDipole(piecen, rslt);
		break;
	case 1: 
		CalcERepulsion(piecen, rslt);
		break;
	case 2: //passed test
		CalcNuclearAttraction(piecen, rslt);
		break;
	case 3: //kinetic
		CalcKinetic(piecen,rslt);
		break;
	case 4: //ions interaction
		CalcIonsInteraction(piecen, rslt);
		break;
	case 5: //passed test
		CalcOverlap(piecen, rslt);
		break;
	default:
		printf("ContractGaussianBasis::GetPrimitiveIntegral() input type=%s is invalide\n",type.c_str());
		printf("Available type are:\n");
		for(std::map<std::string,int>::const_iterator itr=elementType.begin();itr!=elementType.end();++itr){
			printf("\t%s\n",(*itr).first.c_str());
		}
		exit(1);
	}
}

template<class T>
inline T DistanceSqr(const T R1[3],const T R2[3] ){
	return	(R1[0]-R2[0])*(R1[0]-R2[0])+\
			(R1[1]-R2[1])*(R1[1]-R2[1])+\
			(R1[2]-R2[2])*(R1[2]-R2[2]);
}


/*
 Since the recurrence formula does not depend on the function (|r-r'|).
 It is better to implement them in one function.
 */

struct GaussGMFnParameter{
	std::complex<double> T;
	int	GMMaxM;
};

struct GaussFMFnParameter{
	std::complex<double> T;
	int FMMaxM;
};



double GM_auxiFnCos(double u,void * params)
{//auxillary function to calculate GM, which will be called by GSL library.
//calculate Real[u^(2m)*exp(-T*u^2)]
	GaussGMFnParameter * p=static_cast<GaussGMFnParameter *>(params);
	return pow(u*u,p->GMMaxM)*exp(-(p->T.real())*u*u)*cos(p->T.imag()*u*u);
}

double GM_auxiFnSin(double u, void * params)
{//auxillary function to calculate GM, which will be called by GSL library.
//calculate Imag[u^(2m)*exp(-T*u^2)]
	GaussGMFnParameter * p=static_cast<GaussGMFnParameter *>(params);
	return pow(u*u,p->GMMaxM)*exp(-(p->T.real())*u*u)*sin(p->T.imag()*u*u);
}

std::complex<double> ContractGaussianBasis::CalcG_MaxM(std::complex<double> &T){
//get the value of GM(T)
//GM(T)=\sqrt(2/pi)\int_0^1 t^(2m)*exp(-T*t^2)dt;

	const double epsabs=Constants::IntegrationEpsAbs;
	const double epsrel=Constants::IntegrationEpsRel;
	double rslt_r,rslt_i,abserr;
	int level=1000;
	//gsl function
	gsl_function f;
	f.function= &GM_auxiFnCos;
	GaussGMFnParameter gmFnP;
	gmFnP.GMMaxM=GMMaxM;
	gmFnP.T=T;
	f.params=&gmFnP;
	gsl_integration_workspace * workspace=gsl_integration_workspace_alloc (level);
	gsl_integration_qag(&f, 0.0, 1.0, epsabs,epsrel, level, 6,workspace, &rslt_r, &abserr);
	if(abserr>epsabs){
		std::cerr<<"ContractGaussianBasis::CalcGM real part::integration is not converged"<<std::endl;
		gsl_integration_workspace_free(workspace);
		throw "CalcG_MaxM: real part integration is not convereged";
	}
	//imaginary part;
	f.function=&GM_auxiFnSin;
	gsl_integration_qag(&f, 0.0, 1.0, epsabs,epsrel, level, 6,workspace, &rslt_i, &abserr);
	if(abserr>epsabs){
		std::cerr<<"ContractGaussianBasis::CalcGM imag part::integration is not converged"<<std::endl;
		gsl_integration_workspace_free(workspace);
		throw "CalcG_MaxM: imag part integration is not convereged";
	}
	gsl_integration_workspace_free(workspace);
	std::complex<double> rslt=std::complex<double>(rslt_r,-rslt_i);
//	std::complex<double> twoTPowerM=pow(2.0*T,GMMaxM);
	rslt=sqrt(2.0/Constants::PI)*rslt; //GM
	return rslt;
}

bool findPos(const std::list<std::complex<double> > & lst,const std::complex<double> & T, int & idx)
{// look up is not the efficient one.
	bool  rslt=false;
	idx=0;
	for(std::list<std::complex<double> >::const_iterator iter=lst.begin();iter!=lst.end();++iter){
		if(abs((*iter)-T)<1e-12){//T1 exist
			rslt=true;
			break;
		}
		++idx;
	}
	return rslt;
}

inline void ContractGaussianBasis::Insert2GMTable(const std::complex<double> & T, const std::complex<double> & v){
//insert calculated GM(T1)=v into GMVTable and GMTable;
//element is not ordered.
//todo: ordered list might be quicker to find it value.
	GMTTable.push_back(T);
	GMVTable.push_back(v);
}

std::complex<double> ContractGaussianBasis::CalcGM(std::complex<double> & T1,int m)//cause a problem for the calculation of e-repulsion.
{//T=T1^2; T1=\theta*R;
//calculate GM(T)=\sqrt(2/pi)*(2*T)^(m+1/2)*\int_0^1 t^(2m)*exp(-T*t^2)dt;
//G0(\infinity)=1;
//first check if GM(T1) is exist if not, calculate it.
	std::complex<double> rslt;
	double T_critic=1000; //
	std::complex<double> T=T1*T1;
	//bool	epsilon_hm=(-T.real()+m*log(2*(T1.real()*T1.real()+T1.imag()*T1.imag()))-log(1e-12))<0? true:false;//epsilon_hm=true, hm can be omitted
	if(T.real()<T_critic){
		int idx;
		if(findPos(GMTTable,T,idx)){//look up previous calculation
			std::list<std::complex<double> >::iterator iter=GMVTable.begin();
			while(idx>0) {++iter;--idx;}
			rslt=(*iter);
		}else{//G_MaxM not exist
			rslt=CalcG_MaxM(T);
			Insert2GMTable(T,rslt);
		}
		//use recursion relationship to get G_m+1=(2m+1)*G_m-(2T)^m*(sqrt(2)*T1)*exp(-T); not stable
		//use recursion relationship to get G_m=(2*T*G_m+sqrt(2/pi)*exp(-T))/(2*m+1); table
		std::complex<double> h=sqrt(2.0/Constants::PI)*exp(-T);
		for(int i=GMMaxM-1;i>=m;--i){
			rslt=(2.0*T*rslt+h)/(2.0*i+1.0);
		}
	}else{//asymptotic behavior Gm=(2m-1)!!/(2*T)^m*G0=(2m-1)!!/(2*T)^(m+1/2), asymptotic value is not stored.
		rslt=1;
		for(int i=1;i<=m;++i) rslt*=(2.0*i-1.0)/(2.0*T);
		rslt/=sqrt(2.0*T);
	}
	return rslt;
}

std::complex<double> ContractGaussianBasis::CalcM0_Overlap(const QuartetPara & qPara, const int m)
{
	std::complex<double> rslt;
	rslt=Constants::PI/(qPara.xi+qPara.eta)*std::sqrt(Constants::PI/(qPara.xi+qPara.eta))*std::pow(2.0*qPara.thetaSqr,m)*exp(-qPara.thetaSqr*qPara.RSqr);
	return rslt;
}
std::complex<double> ContractGaussianBasis::CalcM0_ERepulsion(const QuartetPara & qPara, const int m)
{
	std::complex<double> rslt;
	std::complex<double> T1=sqrt(qPara.thetaSqr*qPara.RSqr);
	rslt=qPara.GAB*qPara.GCD*std::pow((Constants::PI),3)/(qPara.xi*qPara.eta*std::sqrt(qPara.xi*qPara.eta))*pow(2.0*qPara.thetaSqr,m)*sqrt(2.0*qPara.thetaSqr)*CalcGM(T1,m); //M0
	return rslt;
}

std::complex<double> ContractGaussianBasis::CalcM0_NA(const QuartetPara & qPara, const int m)
{//calculate the integration [0]^m_NA
	std::complex<double> rslt;
	std::complex<double> T1=std::sqrt(qPara.thetaSqr*qPara.RSqr);
	rslt=qPara.GAB*Constants::PI/qPara.xi*std::sqrt(Constants::PI/qPara.xi)*pow(2.0*qPara.thetaSqr,m)*sqrt(2.0*qPara.thetaSqr)*CalcGM(T1,m); //has a problem when R=0;
	return rslt;
}

std::complex<double> ContractGaussianBasis::CalcM0(const QuartetPara & qPara,int m, const std::string  & type)
{//calculate [0]^m
	std::complex<double> rslt;
//form basic shell quartet parameters.
	switch(elementType.at(type))//todo: valide type should be checked
	{
//	case 1://kinetic will be converted to overlap at up level;
	case 5: //overlap
		rslt=CalcM0_Overlap(qPara, m);
		break;
//	case 3://dipole can also be converted to overlap
	case 1://e-repuslion
		rslt=CalcM0_ERepulsion(qPara,m);
		break;
	case 2: //nuclearAttraction
		rslt=CalcM0_NA(qPara,m);
		break;
	default:
		std::cerr<<"Gaussian::CalcM0() invalid input parameter type: "<< type <<std::endl;
		break;
	}
	return rslt;
}

inline void ContractGaussianBasis::GetR(const RM & rm, const int & crntm, int r[3])
{//get the current r based on crntm and the value rm.r;
	int m=crntm;
	if(m<=rm.r[0]){
		r[0]=rm.r[0]-m;
		r[1]=rm.r[1];
		r[2]=rm.r[2];
	}else if(m<=rm.r[0]+rm.r[1]){
		r[0]=0;
		r[1]=rm.r[1]-(m-rm.r[0]);
		r[2]=rm.r[2];
	}else if(m<=rm.r[0]+rm.r[1]+rm.r[2]){
		r[0]=0;
		r[1]=0;
		r[2]=rm.r[2]-(m-rm.r[0]-rm.r[1]);
	}else{
		std::cerr<<"ContractGaussianBasis::GetR:: the value of crntm is not correct"<<std::endl;
	}
}

std::complex<double> ContractGaussianBasis::CalcRM(const QuartetPara &qPara, const RM & rm, const std::string & type)
{
	std::complex<double> rslt;
	std::string t_type;
	switch(elementType.at(type)){
	case 0:
		printf("ContractGaussianBasis::CalcRM() Dipole calculation should be transformed into overlap calculation\n");
		exit(1);
		break;
	case 2:
		rslt=RMKernel(qPara,rm,type);
		break;
	case 1:
		rslt=RMKernel(qPara,rm,type);
		break;
	case 5:
		rslt=RMKernel(qPara,rm,type);
		break;
	case 3://transformed to overlap calculation.
		RM t_rm;
		t_type="overlap";
		memcpy(&t_rm,&rm,sizeof(rm));
		t_rm.r[0]+=2;
		rslt=RMKernel(qPara,t_rm,t_type);

		t_rm.r[0]-=2;t_rm.r[1]+=2;
		rslt+=RMKernel(qPara,t_rm,t_type);

		t_rm.r[1]-=2;t_rm.r[2]+=2;
		rslt+=RMKernel(qPara,t_rm,t_type);

		rslt*=-0.5;
		break;
	default:
		std::cerr<<"ContractGaussianBasis::CalcRM() input type is invalide"<<std::endl;
		break;
	}
	return rslt;
}
std::complex<double> ContractGaussianBasis::RMKernel(const QuartetPara &qPara, const RM & rm, const std::string & type)
{//Calculate the value [r]^m with non-recursive method
	int crntm=rm.r[0]+rm.r[1]+rm.r[2];
	std::vector< std::complex<double> > leafNode(crntm/2+1);//store the value of the current leaves node at line crntm;
	std::complex<double> M0;
	M0=CalcM0(qPara,crntm+rm.m,type);
	leafNode[0]=M0; //only one node for last level.
	--crntm;
	int	r[3];
	while(crntm>=0){
		GetR(rm, crntm,r);//get r value at line crntm-1 of the first left node.
		int ri=0; //indicator from which (x,y,z) begin to deduct.
		while(r[ri]==0) ++ri;
		for(int i=0;i<=crntm && ri<3;i++)
		{//the number of node at line crntm is min{crntm+1, r[0]+r[1]+r[2]}
			switch(r[ri]){
			case	0:
				if(ri==2){//[0]^m
					leafNode[i]=CalcM0(qPara,crntm+rm.m,type);
					++ri;
				}else{//if not [n,0,m]
					++ri;
					--i;//but keep at the origin place.
				}
				break;
			case	1:
				leafNode[i]=qPara.R[ri]*leafNode[i];
				--r[ri];
				ri+=(ri<2)? 1:0;
				break;
			default:
				leafNode[i]=qPara.R[ri]*leafNode[i]-(r[ri]-1.0)*leafNode[i+1];
				--r[ri];
				break;
			}
		}
		--crntm;
	}
	return leafNode[0];
}

template<class T>
inline void VecSub(T R1[3],T R2[3], T R[3]){
	R[0]=R1[0]-R2[0];
	R[1]=R1[1]-R2[1];
	R[2]=R1[2]-R2[2];
}

void ContractGaussianBasis::CalcQuartetPara(const std::vector<PrimitiveGauss * > & ABCD, QuartetPara & qPara ,const std::string & type)
{
	/*
	std::map<std::string, int> elementType;
	elementType["dipole"]=0; // can be converted to overlap
	elementType["eRepulsion"]=1;
	elementType["nuclearAttraction"]=2;
	elementType["kinetic"]=3;//can ve converted to overlap
	elementType["ionsInteraction"]=4;
	elementType["overlap"]=5;
	*/
	std::complex<double> d_R[3];
	std::complex<double> P[3];
	std::complex<double> Q[3];
	switch(elementType.at(type)){
	case 3: //kinetic
	case 5: //overlap
	case 0: //dipole, ABCD[1] and ABCD[3] does not exist
		VecSub(ABCD[2]->R, ABCD[0]->R,d_R);
		qPara.R[0]=d_R[0];
		qPara.R[1]=d_R[1];
		qPara.R[2]=d_R[2];
		qPara.xi=ABCD[0]->alpha;
		qPara.eta=ABCD[2]->alpha;
		qPara.GAB=1.0;
		qPara.GCD=1.0;
		qPara.thetaSqr=qPara.xi*qPara.eta/(qPara.xi+qPara.eta);
		qPara.RSqr=qPara.R[0]*qPara.R[0]+qPara.R[1]*qPara.R[1]+qPara.R[2]*qPara.R[2];
		break;
	case 1: //epulsion
		qPara.xi=ABCD[0]->alpha+ABCD[1]->alpha;
		qPara.eta=ABCD[2]->alpha+ABCD[3]->alpha;
		for(int i=0;i<3;++i){
			P[i]=(ABCD[0]->alpha*ABCD[0]->R[i]+ABCD[1]->alpha*ABCD[1]->R[i])/qPara.xi;
			Q[i]=(ABCD[2]->alpha*ABCD[2]->R[i]+ABCD[3]->alpha*ABCD[3]->R[i])/qPara.eta;
		}
		VecSub(Q,P,qPara.R);
		qPara.GAB=exp(-ABCD[0]->alpha*ABCD[1]->alpha/qPara.xi*DistanceSqr(ABCD[0]->R,ABCD[1]->R));
		qPara.GCD=exp(-ABCD[2]->alpha*ABCD[3]->alpha/qPara.eta*DistanceSqr(ABCD[2]->R,ABCD[3]->R));
		qPara.thetaSqr=qPara.xi*qPara.eta/(qPara.xi+qPara.eta);
		qPara.RSqr=qPara.R[0]*qPara.R[0]+qPara.R[1]*qPara.R[1]+qPara.R[2]*qPara.R[2];
		break;
	case 2: //nuclearAttraction
		qPara.xi=ABCD[0]->alpha+ABCD[1]->alpha; //eta does not exist, and also will not be used;
		//eta is infinity
		for(int i=0;i<3;++i){
			P[i]=(ABCD[0]->alpha*ABCD[0]->R[i]+ABCD[1]->alpha*ABCD[1]->R[i])/qPara.xi;
			Q[i]=ABCD[2]->R[i];
		}
		VecSub(Q,P,qPara.R);
		qPara.GAB=exp(-ABCD[0]->alpha*ABCD[1]->alpha/qPara.xi*DistanceSqr(ABCD[0]->R,ABCD[1]->R));
		qPara.GCD=1.0;
		qPara.thetaSqr=qPara.xi;
		qPara.RSqr=qPara.R[0]*qPara.R[0]+qPara.R[1]*qPara.R[1]+qPara.R[2]*qPara.R[2];
		break;
	default:
		std::cerr<<"ContractGaussianBasis::CalcVecR type value is not valid"<<std::endl;
	}
}

void ContractGaussianBasis::CalcABInt(std::vector<PrimitiveGauss *> AB, std::vector<std::complex<double> >::iterator & beginIter, std::vector<std::complex<double> >::iterator & endIter, const size_t stride)
{//calculate |AB] from |00P]
//the vector is modified, and the returned value is the last one.
	//if only one element there, there is no need to do further calculation
	std::vector<std::complex<double> >::iterator t_beginIter,t_endIter;
	t_beginIter=beginIter;
	advance(t_beginIter,stride);
	if(t_beginIter==endIter)
		return; //only one element, do not need further calculation.

	size_t ab[2][3]={{0}};
	std::complex<double> xi,mCoeff[2][3]={{0}};
	MDTransformInfo mdTransInfo;
	xi=AB[0]->alpha;
	for(size_t i=0;i<3;++i)
		ab[0][i]=AB[0]->L[i];

	if(NULL!=AB[1]){
		xi+=AB[1]->alpha;
		for(size_t i=0;i<3;++i){
			ab[1][i]=AB[1]->L[i];
			mCoeff[0][i]=(AB[1]->R[i]-AB[0]->R[i])/xi;
			mCoeff[1][i]=-mCoeff[0][i]*AB[0]->alpha;
			mCoeff[0][i]*=AB[1]->alpha;
		}
	}

	mdTransInfo.xi=xi;
//form [az,bz]
	int blockPx=(ab[0][1]+ab[1][1]+1)*(ab[0][2]+ab[1][2]+1);
	int blockPy=(ab[0][2]+ab[1][2]+1);
	mdTransInfo.stride=stride;
	for(size_t Px=0; Px<=ab[0][0]+ab[1][0]; ++Px){
		for(size_t Py=0; Py<=ab[0][1]+ab[1][1];++Py){
			t_beginIter=beginIter;
			advance(t_beginIter,(Px*blockPx+Py*blockPy)*stride);
			t_endIter=t_beginIter;
			advance(t_endIter,blockPy*stride);
			//transform to [0Bz| first
			mdTransInfo.mCoeff=mCoeff[1][2];
			mdTransInfo.beginIter=t_beginIter;
			mdTransInfo.endIter=t_endIter;
			mdTransInfo.level=ab[1][2];
			MDTransform1v(mdTransInfo);
			//transform to [AzBz|
			t_beginIter=beginIter;
			advance(t_beginIter,(Px*blockPx+Py*blockPy)*stride);
			t_endIter=t_beginIter;
			advance(t_endIter,blockPy*stride);
			advance(t_beginIter,ab[1][2]*stride);
			mdTransInfo.mCoeff=mCoeff[0][2];
			mdTransInfo.beginIter=t_beginIter;
			mdTransInfo.endIter=t_endIter;
			mdTransInfo.level=ab[0][2];
			MDTransform1v(mdTransInfo);
		}
	}

//form [[ay,az],[by,bz]]
	mdTransInfo.stride=stride*blockPy;
	for(size_t Px=0;Px<=ab[0][0]+ab[1][0];++Px){
		t_beginIter=beginIter;
		advance(t_beginIter,(Px*blockPx+blockPy-1)*stride);
		t_endIter=t_beginIter;
		advance(t_endIter,blockPx*stride);
		//[[0,az],[by,bz]]
		mdTransInfo.mCoeff=mCoeff[1][1];
		mdTransInfo.beginIter=t_beginIter;
		mdTransInfo.endIter=t_endIter;
		mdTransInfo.level=ab[1][1];
		MDTransform1v(mdTransInfo);
		//transform to [[ay,az],[by,bz]]
		t_beginIter=beginIter;
		advance(t_beginIter,(Px*blockPx+blockPy-1)*stride);
		t_endIter=t_beginIter;
		advance(t_endIter,blockPx*stride);//changed
		advance(t_beginIter,ab[1][1]*blockPy*stride);
		mdTransInfo.mCoeff=mCoeff[0][1];
		mdTransInfo.beginIter=t_beginIter;
		mdTransInfo.endIter=t_endIter;
		mdTransInfo.level=ab[0][1];
		MDTransform1v(mdTransInfo);
	}


//form [[ax,ay,az],[bx,by,bz]]
	mdTransInfo.stride=stride*blockPx;
	t_beginIter=beginIter;
	advance(t_beginIter,(blockPx-1)*stride);
	t_endIter=endIter;
	advance(t_endIter,(blockPx-1)*stride);
	//[[0,ay,az],[bx,by,bz]]
	mdTransInfo.mCoeff=mCoeff[1][0];
	mdTransInfo.beginIter=t_beginIter;
	mdTransInfo.endIter=t_endIter;
	mdTransInfo.level=ab[1][0];
	MDTransform1v(mdTransInfo);

	//transform to [[ax,ay,az],[bx,by,bz]]
	t_beginIter=beginIter;
	advance(t_beginIter,(blockPx-1)*stride);
	t_endIter=endIter;
	advance(t_endIter,(blockPx-1)*stride);
	advance(t_beginIter,ab[1][0]*blockPx*stride);
	mdTransInfo.mCoeff=mCoeff[0][0];
	mdTransInfo.beginIter=t_beginIter;
	mdTransInfo.endIter=t_endIter;
	mdTransInfo.level=ab[0][0];
	MDTransform1v(mdTransInfo);
}

void ContractGaussianBasis::OutputParameters()
{//print all related information for check.
//in windows \r has to add as a end of newline.
	if(io->parallel->rank==0){
		std::ofstream * ofile=io->outPut.GetOutputFile("OutputParameters");
		if(NULL==ofile){
			printf("ContractGaussianBasis::OutputParameter() io should set outputparameter file ahead\n");
			exit(1);
		}
		std::string header(15,'*');
		(*ofile)<<header<<"GaussianBasis"<<header<<std::endl;
		*ofile<<"\%GaussianBasis\n";//todo, output the basis set
		char L[]={'S','P','D','F','G','H'};
		for(int n=0;n<basis.size();++n){
			*ofile<<basis[n][0].L[0]<<" "<<basis[n][0].L[1]<<" "<<basis[n][0].L[2]<<"|"
				   <<basis[n][0].R[0]<<" "<<basis[n][0].R[1]<<" "<<basis[n][0].R[2]<<"|";
			for(int c=0;c<basis[n].size();++c){
					*ofile<<basis[n][c].coeff<<"|"<<basis[n][c].alpha<<"|";
			}
			*ofile<<std::endl;
		}	
		*ofile<<"%\n";
		*ofile<<"shellPairCutoff="<<shellPairCutoff<<"\n";
		*ofile<<"signifShellPairRelDist="<<signifShellPairRelDist<<"\n";
		(*ofile)<<std::endl;
	}
}

std::complex<double> ContractGaussianBasis::CalcABCDInt(const std::vector<PrimitiveGauss * > & ABCD, const std::string & type)
{//calculate [AB|CD] with recursion relationship
//ABCD is a four-element vector contains the address of four primitive Gaussian basis
	size_t abcd[4][3]={{0}};
	int Np[4], Nq[4];//tree level for Np and Nq. root node not included.

	for(size_t i=0;i<3;++i){
		abcd[0][i]=ABCD[0]->L[i];
		abcd[2][i]=ABCD[2]->L[i];
	}

	if(NULL!=ABCD[1]){ //for kinetic and overlap, b and d do not exist
		for(int i=0;i<3;++i){
			abcd[1][i]=ABCD[1]->L[i];
		}
	}

	if(NULL!=ABCD[3]){
		for(int i=0;i<3;++i)
			abcd[3][i]=ABCD[3]->L[i];
	}

	for(size_t i=0;i<3;++i){
		Np[i]=abcd[0][i]+abcd[1][i];
		Nq[i]=abcd[2][i]+abcd[3][i];
	}

	Np[3]=(Np[0]+1)*(Np[1]+1)*(Np[2]+1);
	Nq[3]=(Nq[0]+1)*(Nq[1]+1)*(Nq[2]+1);
	//forming a table. Since  [p|q]=(-1)^q[p+q|0]=(-1)^q[r], symmetry should be used to avoid
	QuartetPara qPara;

	CalcQuartetPara(ABCD,qPara,type);
	RM rm;
	rm.m=0;
	std::vector<std::complex<double> > pqTable(Np[3]*Nq[3],std::complex<double>(0,0));
	std::complex<double> rInt;
	int row,col;
	size_t colMax=(Nq[0]+1)*(Nq[1]+1)*(Nq[2]+1);
	size_t rowMax=(Np[0]+1)*(Np[1]+1)*(Np[2]+1);
	int qx,qy,qz;
	for(int rx=0;rx<=Np[0]+Nq[0];++rx) 
		for(int ry=0;ry<=Np[1]+Nq[1];++ry)
			for(int rz=0;rz<=Np[2]+Nq[2];++rz){
				rm.r[0]=rx;rm.r[1]=ry;rm.r[2]=rz;
				rInt=CalcRM(qPara,rm,type); //return [r];
				//find all posible position
				for(int px=0;px<=rx && px<=Np[0];++px)//qx=rx-px;
					for(int py=0;py<=ry && py<=Np[1];++py)//qy=ry-py
						for(int pz=0;pz<=rz && pz<=Np[2];++pz){//dimension is([px][py][qz])x([qx][qy][qz)
							qx=rx-px;
							qy=ry-py;
							qz=rz-pz;
							if(qx<=Nq[0] && qy<=Nq[1] && qz<=Nq[2]){
								row=(px*(Np[1]+1)+py)*(Np[2]+1)+pz;
								col=(qx*(Nq[1]+1)+qy)*(Nq[2]+1)+qz;
								pqTable[row*colMax+col]=rInt*pow(-1,qx+qy+qz);
							}
						}
			}


	//calculate CD first;
	std::vector<std::complex<double> >::iterator t_beginIter,t_endIter,beginIter,endIter;
	beginIter=pqTable.begin();
	endIter=pqTable.end();
	std::vector<PrimitiveGauss *> AB(2,static_cast<PrimitiveGauss *>(NULL));

	if(colMax>1){//if only one column exist, we do not need further calculation.
		AB[0]=ABCD[2]; AB[1]=ABCD[3];
		t_beginIter=beginIter;
		for(size_t rn=0;rn<rowMax;++rn){
			t_endIter=t_beginIter;
			advance(t_endIter,colMax);
			CalcABInt(AB,t_beginIter,t_endIter,1);
			t_beginIter=beginIter;
			advance(t_beginIter,colMax*(rn+1));
		}
	}

	//use the last column to calculate AB
	if(rowMax>1){//if only one row, we also do not need further calculation.
		t_beginIter=beginIter;
		advance(t_beginIter,colMax-1);
		t_endIter=endIter;
		advance(t_endIter,colMax-1);
		AB[0]=ABCD[0];
		AB[1]=ABCD[1];
		CalcABInt(AB,t_beginIter,t_endIter,colMax);
	}
	return pqTable.back();
}

void ContractGaussianBasis::MDTransform2v(const MDTransformInfo & mdTransInfo)
{//Do MD transform for the array from beiginIter to endIter with stride.
//1v means, only oone variable each time transform from [abp], and p always start
//from zero and increase linearly in the vector. At the end, we have [a,b+p];
//for example, if we want to get [ab], our input vector should be
//[a00],[a01],[a02],...[a0b], and the output will be
//[a0],[a1],[a2],...[ab].
//while 2v means, we need get all possible integral for [axb]. So
//p does not increase linearly.
//For example, if we need to get all possible [axb]. input vector should be
//[000][001][002]...[00b][001][002]...[00(b+1)][002][003]...[00(b+2)][003][004]....
//the output will be
//[00][01][02]...[0b][10][11]..[1b][20][21]...[2b][30][31]....
//the value a,b is set in the variable mdTransInfo.ab[2];
//MDTransformInfo.level is not useful in such case.

	std::vector<std::complex<double> >::iterator t_beginIter,t_endIter,iter;
	//check if there is only one element there
	t_beginIter=mdTransInfo.beginIter;
	advance(t_beginIter,mdTransInfo.stride);
	if(t_beginIter==mdTransInfo.endIter){ //only one element exist.
		return;
	}

//iteration for a
	std::complex<double> t_eles[3];
	std::complex<double> Iover2xi=1.0/(2.0*mdTransInfo.xi);
	std::vector<std::complex<double> > t_vec(mdTransInfo.ab[0]+mdTransInfo.ab[1]+1,0);

	//copy element into t_vec;
	t_beginIter=mdTransInfo.beginIter;
//	t_endIter=t_beginIter;
//	advance(t_endIter,mdTransInfo.stride*(mdTransInfo.ab[1]+1));
	size_t n=0;
	while(n<mdTransInfo.ab[0]){
		t_vec[n]=*(t_beginIter);
		n++;
		advance(t_beginIter,mdTransInfo.stride*(mdTransInfo.ab[1]+1));
	}
	//copy [00(a)],...[00(a+b)]
	n=0;
	while(n<=mdTransInfo.ab[1]){
		t_vec[mdTransInfo.ab[0]+n]=*t_beginIter;
		advance(t_beginIter,mdTransInfo.stride);
		++n;
	}

	//iteration to calculate a
	double p;
	t_beginIter=mdTransInfo.beginIter;
	advance(t_beginIter,mdTransInfo.stride*(mdTransInfo.ab[1]+1));
	for(size_t n=1;n<=mdTransInfo.ab[0];++n){//after each iteration, a+1
		t_eles[0]=0;
		t_eles[1]=t_vec[n-1];
		p=0;
		for(size_t m=n;m<t_vec.size();++m){
			t_eles[2]=t_vec[m];
			t_vec[m]=p*t_eles[0]+mdTransInfo.mCoeffa*t_eles[1]+Iover2xi*t_eles[2];
			t_eles[0]=t_eles[1];
			t_eles[1]=t_eles[2];
			p+=1.0;
		}
		//store [n00],[n01],...,[n0b] to the value to nth segment
		for(size_t m=0;m<=mdTransInfo.ab[1];++m){
			*t_beginIter=t_vec[m+n];
			advance(t_beginIter,mdTransInfo.stride);
		}
	}

	//for each segmentation of a, calculate b
	MDTransformInfo t_mdTransInfo;
	if(mdTransInfo.ab[1]>0){
		t_mdTransInfo.xi=mdTransInfo.xi;
		t_mdTransInfo.stride=mdTransInfo.stride;
		t_mdTransInfo.level=mdTransInfo.ab[1];
		t_mdTransInfo.mCoeff=mdTransInfo.mCoeff; //coeffecient for b;
		t_mdTransInfo.beginIter=mdTransInfo.beginIter;
		for(size_t m=0;m<=mdTransInfo.ab[0];++m){
			t_mdTransInfo.endIter=t_mdTransInfo.beginIter;
			advance(t_mdTransInfo.endIter,mdTransInfo.stride*(mdTransInfo.ab[1]+1));
			MDTransform1v(t_mdTransInfo);
			t_mdTransInfo.beginIter=t_mdTransInfo.endIter;
		}
	}

}

void ContractGaussianBasis::MDTransform1v(const MDTransformInfo & mdTransInfo)
{//Do MD transform for the array from beiginIter to endIter with stride.
//1v means, only oone variable each time transform from [abp], and p always start
//from zero and increase linearly in the vector. At the end, we have [a,b+p];
//for example, if we want to get [ab], our input vector should be
//[a00],[a01],[a02],...[a0b], and the output will be
//[a0],[a1],[a2],...[ab].
//while 2v means, we need get all possible integral for [axb]. So
//p does not increase linearly.
//For example, if we need to get all possible [axb]. input vector should be
//[000][001][002]...[00b][001][002]...[00(b+1)][002][003]...[00(b+2)][003][004]....
//the output will be
//[00][01][02]...[0b][10][11]..[1b][20][21]...[2b][30][31]....
	std::vector<std::complex<double> >::iterator iter,t_beginIter,t_endIter;
	//check if only one elelment exist
	t_beginIter=mdTransInfo.beginIter;
	advance(t_beginIter,mdTransInfo.stride);
	if(t_beginIter==mdTransInfo.endIter){
		return;
	}


	std::complex<double> t_eles[3];
	std::complex<double> Iover2xi=1.0/(2.0*mdTransInfo.xi);
	double p=0;
	t_beginIter=mdTransInfo.beginIter;
	t_endIter=t_beginIter;
	advance(t_endIter,mdTransInfo.level*mdTransInfo.stride);//when to stop the iteration
	while(t_beginIter!=t_endIter){
		p=-1;
		t_eles[0]=0;
		t_eles[1]=*t_beginIter;
		advance(t_beginIter,mdTransInfo.stride);
		iter=t_beginIter;
		while(iter!=mdTransInfo.endIter){
			p+=1.0;
			t_eles[2]=*iter;// can only put here, since we need to check whether iter is valid.
			*iter=p*t_eles[0]+mdTransInfo.mCoeff*t_eles[1]+Iover2xi*t_eles[2];
			advance(iter,mdTransInfo.stride);
			t_eles[0]=t_eles[1];
			t_eles[1]=t_eles[2];
		}
	}
}

void ContractGaussianBasis::MDTransformVector(const std::vector<PrimitiveGauss *> & CD,  std::vector<std::complex<double> >::iterator & beginIter,std::vector<std::complex<double> >::iterator & endIter, const size_t stride)
{//calculate |CD] with recursion relationship
//begin and end iterator point to the initial value
//calculated value is also stored in the vector to replace original data
//data organized as [cx][dx][cy][dy][cz][dz]
	std::vector<std::complex<double> >::iterator t_beginIter,t_endIter;

	t_beginIter=beginIter;
	advance(t_beginIter,stride);
	if(t_beginIter==endIter){
		return;
	}

	size_t cd[2][3]={{0}};
	std::complex<double> xi, mCoeff[2][3]={{std::complex<double>(0,0)}};//std:complex<double>(0,0);
	xi=CD[0]->alpha;

	for(int i=0;i<3;++i)
		cd[0][i]=CD[0]->L[i];

	if(NULL!=CD[1]){//if d exist
		xi+=CD[1]->alpha;
		for(int i=0;i<3;++i){
			mCoeff[0][i]=(CD[1]->R[i]-CD[0]->R[i])/xi;
			mCoeff[1][i]=-mCoeff[0][i]*CD[0]->alpha;
			mCoeff[0][i]*=CD[1]->alpha;
			cd[1][i]=CD[1]->L[i];
		}
	}

	MDTransformInfo mdTransInfo;
	mdTransInfo.xi=xi; //is the same for all dimensions.
	//recursion calculate cx and dx at each block. while each block has cy*dy*cz*dz elements.
	//chose nth element from each block, to form a new array and then do recursion process.
	//after whole iteration, we get [[cx,dx]]
	size_t blockSizeCxDx=(cd[0][1]+1)*(cd[0][2]+1)*(cd[1][1]+1)*(cd[1][2]+1);
	if(cd[0][0]+cd[1][0]>0){//if it is zero, there is only one element exist, no need to do further calculation
		mdTransInfo.ab[0]=cd[0][0];
		mdTransInfo.ab[1]=cd[1][0];
		mdTransInfo.stride=stride*blockSizeCxDx;
		mdTransInfo.mCoeff=mCoeff[1][0];
		mdTransInfo.mCoeffa=mCoeff[0][0];
		for(size_t i=0;i<blockSizeCxDx;++i){
			t_beginIter=beginIter;
			advance(t_beginIter,i*stride);
			t_endIter=endIter;
			advance(t_endIter,i*stride);
			mdTransInfo.beginIter=t_beginIter;
			mdTransInfo.endIter=t_endIter;
			MDTransform2v(mdTransInfo);
		}
	}

	//recursion to calculate cy, dy  at each block and each block has a size of cz*dz;
	size_t blockSizeCyDy=(cd[0][2]+1)*(cd[1][2]+1);
	if(cd[0][1]+cd[1][1]>0){
		mdTransInfo.stride=blockSizeCyDy*stride;
		mdTransInfo.ab[0]=cd[0][1];
		mdTransInfo.ab[1]=cd[1][1];
		mdTransInfo.mCoeffa=mCoeff[0][1];
		mdTransInfo.mCoeff=mCoeff[1][1];
		for(size_t cx=0;cx<=cd[0][0];++cx)
			for(size_t dx=0;dx<=cd[1][0];++dx)
				for(size_t i=0;i<blockSizeCyDy;++i){
					t_beginIter=beginIter;
					advance(t_beginIter,((cx*(cd[1][0]+1)+dx)*blockSizeCxDx+i)*stride);
					t_endIter=t_beginIter;
					advance(t_endIter,blockSizeCxDx*stride);
					mdTransInfo.beginIter=t_beginIter;
					mdTransInfo.endIter=t_endIter;
					MDTransform2v(mdTransInfo);
				}
	}
	//recursion to calculate cz,dz at each block and each block has a size of ;
	if(cd[0][2]+cd[1][2]>0){
		mdTransInfo.stride=stride;
		mdTransInfo.ab[0]=cd[0][2];
		mdTransInfo.ab[1]=cd[1][2];
		mdTransInfo.mCoeff=mCoeff[1][2];
		mdTransInfo.mCoeffa=mCoeff[0][2];
		for(size_t cx=0;cx<=cd[0][0];++cx)
			for(size_t dx=0;dx<=cd[1][0];++dx)
				for(size_t cy=0;cy<=cd[0][1];++cy)
					for(size_t dy=0;dy<=cd[1][1];++dy){
							t_beginIter=beginIter;
							advance(t_beginIter,((cx*(cd[1][0]+1)+dx)*blockSizeCxDx+(cy*(cd[1][1]+1)+dy)*blockSizeCyDy)*stride);
							t_endIter=t_beginIter;
							advance(t_endIter,blockSizeCyDy*stride);
							mdTransInfo.beginIter=t_beginIter;
							mdTransInfo.endIter=t_endIter;
							MDTransform2v(mdTransInfo);
						}
	}
}

void ContractGaussianBasis::MDTransformTable(const std::vector<PrimitiveGauss *> & ABCD, std::vector<std::complex<double> > & pqTable)
{//calcualte [AB|CD] from [p|q]
//the output table contains all integrations [ab|cd] where a={0,1,..,A}, b={0,1,...,B},c={0,1,...,C},d={0,1,...,D}
//idxr and idxc are the fundamental gaussian parameters for [AB| and |CD] separately
//pqTable is stored in one dimension
	size_t Pmax=1;
	size_t Qmax=1;
	std::vector<PrimitiveGauss *> AB;
	std::vector<PrimitiveGauss *> CD;
	AB.push_back(ABCD[0]);
	AB.push_back(ABCD[1]);
	CD.push_back(ABCD[2]);
	CD.push_back(ABCD[3]);

	for(int i=0;i<3;++i){
		Pmax*=(AB[0]->L[i]+1);
		Qmax*=(CD[0]->L[i]+1);
	}

	if(NULL!=AB[1]){//if b exist
		for(int i=0;i<3;++i)
			Pmax*=(AB[1]->L[i]+1);
	}

	if(NULL!=CD[1]){// if d exist
		for(int i=0;i<3;++i)
			Qmax*=(CD[1]->L[i]+1);
	}

	//MDTransform for each row.
	std::vector<std::complex<double> >::iterator beginIter, endIter;
	for(size_t i=0;i<Pmax;++i){
		beginIter=pqTable.begin();
		advance(beginIter, i*Qmax);
		endIter=beginIter;
		advance(endIter, Qmax);
		MDTransformVector(CD,beginIter,endIter,1);
	}

	//MDTransform for each column
	for(size_t i=0;i<Qmax;++i){
		beginIter=pqTable.begin();
		advance(beginIter,i);
		endIter=pqTable.end();
		advance(endIter, i);
		MDTransformVector(AB, beginIter, endIter, Qmax);
	}


}

void ContractGaussianBasis::GetFullPQTable(const std::vector<PrimitiveGauss *> &ABCD, std::vector<std::complex<double> > & pqTable,const std::string & type)
{//get a full pqTable which can be used to calcualte the
//integration of [ab|cd] for each a,b, c,and d.
//pqTable will be resize
	int abcd[4][3]={{0}};
	for(size_t i=0;i<3;++i){
		abcd[0][i]=ABCD[0]->L[i];
		abcd[2][i]=ABCD[2]->L[i];
	}

	if(NULL!=ABCD[1]){
		for(size_t i=0;i<3;++i)
			abcd[1][i]=ABCD[1]->L[i];
	}
	if(NULL!=ABCD[3]){
		for(size_t i=0;i<3;++i)
			abcd[3][i]=ABCD[3]->L[i];
	}

	RM rm;
	std::complex<double> rInt;
	QuartetPara	qPara;
	CalcQuartetPara(ABCD,qPara,type);
	rm.m=0;

	size_t rxMax,ryMax,rzMax;
	rxMax=abcd[0][0]+abcd[1][0]+abcd[2][0]+abcd[3][0];
	ryMax=abcd[0][1]+abcd[1][1]+abcd[2][1]+abcd[3][1];
	rzMax=abcd[0][2]+abcd[1][2]+abcd[2][2]+abcd[3][2];
	std::complex<double> r[rxMax+1][ryMax+1][rzMax+1];

	for(size_t rx=0;rx<=rxMax;++rx)
		for(size_t ry=0;ry<=ryMax;++ry)
			for(size_t rz=0;rz<=rzMax;++rz){
				rm.r[0]=rx;
				rm.r[1]=ry;
				rm.r[2]=rz;
				r[rx][ry][rz]=CalcRM(qPara,rm,type);
//				std::cout<<"r["<<rx<<"]["<<ry<<"]["<<rz<<"]="<<r[rx][ry][rz]<<std::endl;
			}

	//iterate around the table and fill it with appropriate [r] value
	//the index of pqtable is [ax][bx][ay][by][az][bz][cx][dx][cy][dy][cz][dz]
	//todo: need to check the value in blocksize
	size_t blocksize[13]={0};
	int	n=12;
	blocksize[12]=1;
	for(int lr=1;lr>=0;--lr) //left(AB) or right(CD)
		for(int x=2;x>=0;--x) //x,y,z
			for(int i=1;i>=0;--i)// which one of AB or CD
			{
				--n;
				blocksize[n]=blocksize[n+1]*(abcd[lr*2+i][x]+1);
			}

	size_t idx[12];
	pqTable.resize(blocksize[0],std::complex<double>(0,0));//allocate enough memory
	int p[3],q[3];
	for(size_t pos=0,t_pos=0;pos<pqTable.size();++pos){
		t_pos=pos;
		for(int i=1;i<=12;++i){
			if(t_pos>=blocksize[i]){
				idx[i-1]=t_pos/blocksize[i];
				t_pos-=idx[i-1]*blocksize[i];
			}else{
				idx[i-1]=0;
			}
//			std::cout<<"idx["<<i-1<<"]="<<idx[i-1]<<std::endl;
		}
		//calculate corresponding p and q
		p[0]=idx[0]+idx[1];
		p[1]=idx[2]+idx[3];
		p[2]=idx[4]+idx[5];
		q[0]=idx[6]+idx[7];
		q[1]=idx[8]+idx[9];
		q[2]=idx[10]+idx[11];
//		std::cout<<"pqTable:"<<"r["<<p[0]+q[0]<<"]["<<p[1]+q[1]<<"]["<<p[2]+q[2]<<"]="<<r[p[0]+q[0]][p[1]+q[1]][p[2]+q[2]]<<std::endl;
		pqTable[pos]=r[p[0]+q[0]][p[1]+q[1]][p[2]+q[2]]*pow(-1,q[0]+q[1]+q[2]);
	}

}

int ContractGaussianBasis::CalcKinetic(int piecen, std::vector<std::complex<double> > & rslt)
{//calculate the kinetic matrix  term;
	int pos_i=piecen*segmentSize;
	int pos_f=pos_i+segmentSize;	
	if(pos_f>localDataSize){
		pos_f=localDataSize;
		rslt.resize((pos_f-pos_i),0);		
	}else{
		rslt.resize(segmentSize,0);		
	}
	int row,col;
	int n=0;	
	std::vector<PrimitiveGauss *> abcd(4,NULL);

	std::complex<double> tmp;
	int rank(io->parallel->rank);
	for(;pos_i<pos_f;++pos_i){
		row=pos_i/matrixSize;
		col=pos_i-row*matrixSize;
		std::vector<PrimitiveGauss> & cgb_r(basis[row+basisIndexRange[rank]]);
		std::vector<PrimitiveGauss> & cgb_c(basis[col]);
		rslt[n]=0;
		for(int r=0;r<cgb_r.size();++r){
			for(int c=0;c<cgb_c.size();++c){
				abcd[0]=&cgb_r[r];	
				abcd[2]=&cgb_c[c];
				tmp=CalcABCDInt(abcd,"kinetic");
				rslt[n]+=tmp*abcd[0]->coeff*abcd[2]->coeff*abcd[0]->norm*abcd[2]->norm;
			}
		}	
		rslt[n]*=rotExpAlpha*rotExpAlpha;
		++n;
	}

	return 0;
}


int ContractGaussianBasis::CalcOverlap(int piecen, std::vector<std::complex<double> > & rslt)
{
	int pos_i=piecen*segmentSize;
	int pos_f=pos_i+segmentSize;	
	if(pos_f>localDataSize){
		pos_f=localDataSize;
		rslt.resize((pos_f-pos_i),0);		
	}else{
		rslt.resize(segmentSize,0);		
	}
	int row,col;
	int n=0;	
	std::vector<PrimitiveGauss *> abcd(4,NULL);
	int rank(io->parallel->rank);
	for(;pos_i<pos_f;++pos_i){
		row=pos_i/matrixSize;
		col=pos_i-row*matrixSize;
		std::vector<PrimitiveGauss> & cgb_r(basis[row+basisIndexRange[rank]]);
		std::vector<PrimitiveGauss> & cgb_c(basis[col]);
		rslt[n]=0;
		for(int r=0;r<cgb_r.size();++r){
			for(int c=0;c<cgb_c.size();++c){
				abcd[0]=&cgb_r[r];	
				abcd[2]=&cgb_c[c];
				rslt[n]+=CalcABCDInt(abcd, "overlap")*abcd[0]->coeff*abcd[2]->coeff*abcd[0]->norm*abcd[2]->norm;
			}
		}	
		++n;
	}
	return 0;
}

int ContractGaussianBasis::CalcIonsInteraction(int piecen, std::vector<std::complex<double> > & rslt){
	CalcOverlap(piecen, rslt);	
	for(int n=0;n<rslt.size();++n){
		rslt[n]*=ionsInteractionEnergy;
	}	
	return 0;
}



int ContractGaussianBasis::CalcNuclearAttraction(int piecen, std::vector<std::complex<double> >& rslt)
{
	int pos_i=piecen*segmentSize;
	int pos_f=pos_i+segmentSize;	
	if(pos_f>localDataSize){
		pos_f=localDataSize;
		rslt.resize((pos_f-pos_i),0);		
	}else{
		rslt.resize(segmentSize,0);		
	}
	int row,col;
	int n=0;	
	std::vector<PrimitiveGauss *> abcd(4,NULL);
	PrimitiveGauss	pgo_c;
	pgo_c.alpha=0;
	pgo_c.L[0]=0;
	pgo_c.L[1]=0;
	pgo_c.L[2]=0;
	abcd[2]=&pgo_c;
	double eleNum;
	std::complex<double> temp;
	int	rank(io->parallel->rank);
	for(;pos_i<pos_f;++pos_i){
		row=pos_i/matrixSize;
		col=pos_i-row*matrixSize;
		std::vector<PrimitiveGauss> & cgb_r(basis[row+basisIndexRange[rank]]);
		std::vector<PrimitiveGauss> & cgb_c(basis[col]);
		rslt[n]=0;
		for(int r=0;r<cgb_r.size();++r){
			for(int c=0;c<cgb_c.size();++c){
				abcd[0]=&cgb_r[r];	
				abcd[1]=&cgb_c[c];
				for(int ni=0;ni<atomLoc->size();++ni){
					pgo_c.R[0]=(*atomLoc)[ni].R[0];
					pgo_c.R[1]=(*atomLoc)[ni].R[1];
					pgo_c.R[2]=(*atomLoc)[ni].R[2];
					eleNum=-(species->species[species->name2idx[(*atomLoc)[ni].name]]).eleNum;//minus comes from electron
					rslt[n]+=CalcABCDInt(abcd, "nuclearAttraction")*abcd[0]->coeff*abcd[1]->coeff*eleNum*abcd[0]->norm*abcd[1]->norm;

				}
			}
		}
		rslt[n]*=rotExpAlpha;	
		++n;
	}
	return 0;
}

int ContractGaussianBasis::GetERepulsionPieceNum(int & dataSegmentSize){
	int eRepulsionLocalDataSize=localDataSize*size*size;
	if(eRepulsionLocalDataSize>dataSegmentSize){
		ee_segmentSize=dataSegmentSize;
		return ceil(eRepulsionLocalDataSize*1.0/dataSegmentSize);
	}else{
		ee_segmentSize=eRepulsionLocalDataSize;
		return 1;
	}
	return 0;
}


int ContractGaussianBasis::CalcERepulsion(int piecen, std::vector<std::complex<double> > & rslt)
{//to be determined later by how to partition the calcualtion
	int pos_i=piecen*ee_segmentSize;
	int pos_f=pos_i+ee_segmentSize;	
	if(pos_f>localDataSize*size*size){//local datasize is not currect since one electron and two electron system, this is different
		pos_f=localDataSize*size*size;
		rslt.resize((pos_f-pos_i),0);		
	}else{
		rslt.resize(ee_segmentSize,0);		
	}
	int row,col,row1,row2,col1,col2;
	int n=0;	
	int eMatrixSize=matrixSize*size;
	std::vector<PrimitiveGauss *> abcd(4,NULL);
	int rank(io->parallel->rank);
	for(;pos_i<pos_f;++pos_i){
		row=pos_i/(eMatrixSize);
		col=pos_i-row*eMatrixSize;
		row1=row/size;
		row2=row-row1*size;
		col1=col/size;
		col2=col-col1*size;	
		std::vector<PrimitiveGauss> & cgb_r1(basis[row1+basisIndexRange[rank]]);
		std::vector<PrimitiveGauss> & cgb_r2(basis[row2]);
		std::vector<PrimitiveGauss> & cgb_c1(basis[col1]);
		std::vector<PrimitiveGauss> & cgb_c2(basis[col2]);
		rslt[n]=0;
		for(int r1=0;r1<cgb_r1.size();++r1){
			for(int r2=0;r2<cgb_r2.size();++r2){
				for(int c1=0;c1<cgb_c1.size();++c1){
					for(int c2=0;c2<cgb_c2.size();++c2){
						abcd[0]=&cgb_r1[r1];	
						abcd[1]=&cgb_r2[r2];	
						abcd[2]=&cgb_c1[c1];
						abcd[3]=&cgb_c2[c2];
						rslt[n]+=CalcABCDInt(abcd, "eRepulsion")*abcd[0]->coeff*abcd[2]->coeff*abcd[0]->norm*abcd[2]->norm;
					}
				}
			}
		}	
		++n;
	}
	return 0;
}


int ContractGaussianBasis::CalcDipole(int piecen, std::vector<std::complex<double> > & rslt)
{//The calculation of <|r|> is transformed to the calcualtion of overlap matrix elements
//[r]=[ab|(c+1)d]+C[ab|cd]
	int pos_i=piecen*segmentSize;
	int pos_f=pos_i+segmentSize;	
	if(pos_f>localDataSize){
		pos_f=localDataSize;
		rslt.resize((pos_f-pos_i)*3,0);		
	}else{
		rslt.resize(segmentSize*3,0);		
	}
	int row,col;
	int n=0;	
	std::complex<double> overlap;
	std::complex<double> overlap_1;
	std::vector<PrimitiveGauss *> abcd(4,NULL);
	std::complex<double> expAlpha=1.0/rotExpAlpha;
	int rank(io->parallel->rank);
	for(;pos_i<pos_f;++pos_i){
		row=pos_i/matrixSize;
		col=pos_i-row*matrixSize;
		std::vector<PrimitiveGauss> & cgb_r(basis[row+basisIndexRange[rank]]);
		std::vector<PrimitiveGauss>         cgb_c(basis[col]);//because we have to change the momentum of cgb_c, so we cannot use reference any more
		rslt[3*n]=0; 
		rslt[3*n+1]=0;
		rslt[3*n+2]=0;
		for(int r=0;r<cgb_r.size();++r){
			for(int c=0;c<cgb_c.size();++c){
				abcd[0]=&cgb_r[r];	
				abcd[2]=&cgb_c[c];
				overlap=CalcABCDInt(abcd,"overlap");
				for(int i=0;i<3;++i){
					++(abcd[2]->L[i]);	
					overlap_1=CalcABCDInt(abcd,"overlap");
					--(abcd[2]->L[i]);
					rslt[3*n+i]+=(overlap_1+overlap*abcd[2]->R[i])*abcd[0]->coeff*abcd[2]->coeff*abcd[0]->norm*abcd[2]->norm;
				}
			}
		}	
		rslt[3*n]*=expAlpha;
		rslt[3*n+1]*=expAlpha;
		rslt[3*n+2]*=expAlpha;
		++n;
	}
	return 0;
}


int  ContractGaussianBasis::GetWaveFunction(const std::vector<double> & xyz,const std::vector< std::complex <double> >::iterator coeffBegin,
					std::vector< std::complex<double> >::iterator wfsBegin){

	return 0;
}

ContractGaussianBasis::~ContractGaussianBasis(){
}


void ContractGaussianBasis::ParseGaussianBasis(){
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

void ContractGaussianBasis::AddNewBasis(const std::vector<int> & L,const std::vector<std::complex<double> > & alpha, const std::vector<std::vector<double> > & coeff, double R[3]){
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


void ContractGaussianBasis::AdjustBasisSize(const std::vector<int> & L){
	int num=0;
	for(int n=0;n<L.size();++n){
			num+=(L[n]==0)? 1:3*L[n];
	}
	basis.resize(basis.size()+num);
}

void ContractGaussianBasis::CalcPrimitiveGaussNorm(){
//after reading Gaussian basis set, calc the norm of each primitive gaussian basis
	int * L;	
	for(int n=0;n<basis.size();++n){
		for(int m=0;m<basis[n].size();++m){
			L=basis[n][m].L;
			basis[n][m].norm=(std::pow(2.0*basis[n][m].alpha,(2.0*(L[0]+L[1]+L[2])+3.0)/4.0))/sqrt(tgamma(L[0]+0.5)*tgamma(L[1]+0.5)*tgamma(L[2]+0.5));			
		}
	}
}

void ContractGaussianBasis::GetBasisSpreadRange(int basis_id,double epsilon, double center[3], double & radii){
//epsilon is the precision
//the function evaluate the spread range of the gaussian basis
	double r(0);
	double log_epsilon(log(epsilon));
	ContractGauss & pb(basis[basis_id]);
	radii=0.0;
	for(int n=0;n<basis[basis_id].size();++n){
		r=(-log_epsilon+0.5*log(pb[n].alpha.real()))/pb[n].alpha.real();
		if(r>radii){
			radii=r;
		}
	}
	radii=sqrt(radii);
	center[0]=pb[0].R[0].real();
	center[1]=pb[0].R[1].real();
	center[2]=pb[0].R[2].real();
}

int  ContractGaussianBasis::CalcBasisValueAtGrids(double xyz[3], std::list<int> basis_id, double * basis_value){
	ContractGauss * p_CGB;	
	std::complex<double> val;
	std::complex<double> *R;
	int * L;
	for(std::list<int>::iterator iter=basis_id.begin();iter!=basis_id.end();++iter){
		p_CGB=&basis[*iter];	
		val=0;
		for(int n=0;n<(*p_CGB).size();++n){
			R=(*p_CGB)[n].R;
			L=(*p_CGB)[n].L;
			val+=(*p_CGB)[n].norm*(*p_CGB)[n].coeff*\
				std::pow(xyz[0]-R[0],L[0])*std::pow(xyz[1]-R[1],L[1])*std::pow(xyz[2]-R[2],L[2])*\
				exp(-((xyz[0]-R[0])*(xyz[0]-R[0])+(xyz[1]-R[1])*(xyz[1]-R[1])+(xyz[2]-R[2])*\
					(xyz[2]-R[2])));
		}
	}
}

#endif
