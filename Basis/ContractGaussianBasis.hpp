#ifndef _GAUSSIANBASIS_H_
#define _GAUSSIANBASIS_H_

#include<vector>
#include"../Constants/Constants.h"
#include<complex>
#include"../Common/Header.h"

class ContractGaussianBasis{

//structure
public:
	struct RM{
		int r[3];
		int m;
	};

	typedef struct{ //Primitive Gaussian Basis.
		std::complex<double> coeff;
		std::complex<double> alpha;
		std::complex<double> R[3];
		int					 L[3];
		std::complex<double> norm;
	} PrimitiveGauss;
	
	typedef std::vector<PrimitiveGauss> ContractGauss;
	
	struct QuartetPara{
		std::complex<double> xi;
		std::complex<double> eta;
		std::complex<double> R[3];
		std::complex<double> GAB, GCD;
		std::complex<double> thetaSqr;
		std::complex<double> RSqr;
	};

	struct MDTransformInfo{
		//with a given array start with beginIter and end with endIter. calculation start from the second element.
		//the start p may not necessarily be zero. if not zero, then vp_2 has to be provided for
		std::vector<std::complex<double> >::iterator beginIter, endIter;
		size_t 				 stride; //stride of the vector.
		std::complex<double> mCoeff;
		std::complex<double> xi;
		size_t 				 level; //iteration of recursion should be done
		size_t				 ab[2];//for function MDTransform2v, ab[0]=a,ab[1]=b;
		std::complex<double> mCoeffa; //for function MDTransform2v, for ab[0]
	};
//data
public:
	//std::map<PrimitiveGauss, int > g2sIdx; //a map from GaussianBasis index to space basis index.
//	std::vector<PrimitiveGauss>		s2gIdx; //a map from general index to Gaussian basis.
//	std::vector<PrimitiveGauss> 	par; //exponential coefficient in Gaussian basis
	double							shellPairCutoff;// not used yet.
	double			  				signifShellPairRelDist;
	std::vector<std::vector<int> >  shellPair;

	std::vector<ContractGauss>		basis;

private:
	int 								GMMaxM;
	std::list< std::complex<double> > 	GMVTable; //for GM calculation
	std::list< std::complex<double> >	GMTTable;//since T is complex, interpolation is complicated, not considered by now.
	std::vector<int>					accummulatedElement;
	
//method
public:
	void help();
	void					GetBasisSpreadRange(int basis_id, double epsilon, double center[3], double & radii);
	int 					Init(IO * io,Species * species, size_t dim);
	bool Check();
	void OutputParameters(); //input parser will output parameters once read.
	int  GetWaveFunction(const std::vector<double> & xyz,
						const std::vector< std::complex <double> >::iterator coeffBegin,
						std::vector< std::complex<double> >::iterator wfsBegin); //get the wfs's valueh at point xyz. The wfs's coefficient is stored in coeff.
	void GetLocalSpaceTermPiecen(int piecen,  const std::string & type,std::vector<std::complex<double> > & rslt);
	int  GetERepulsionPieceNum(int & dataSegmentSize);
	int  CalcBasisValueAtGrids(double xyz[3], std::list<int> basis_id, double * basis_value);


							 ContractGaussianBasis();
							~ContractGaussianBasis();
private:
	bool CheckVar(const std::string &);
	void GetR(const RM & rm, const int & crntm,int r[3]);
	void InputParameters(); //load parameters from input parser.
	int	 CalcKinetic(int piecen, std::vector<std::complex<double> > & rslt);
	int  CalcNuclearAttraction(int piecen, std::vector<std::complex<double> > & rslt);
	int	 CalcDipole(int piecen, std::vector<std::complex<double> > & rslt);
	int  CalcOverlap(int piecen, std::vector<std::complex<double> > & rslt);
	int  CalcIonsInteraction(int piecen, std::vector<std::complex<double> > & rslt);
private: //should be private
	void 				 Insert2GMTable(const std::complex<double> & T1, const std::complex<double> & v);
	std::complex<double> CalcM0(const QuartetPara & qPara, int m, const std::string & type);
	std::complex<double> CalcG_MaxM(std::complex<double> &T1);
	std::complex<double> CalcGM(std::complex<double> & T,int m);//GM definition is \sqrt(2/pi)*(2*T)^(m+1/2)*\int_0^1 t^(2m)*exp(-T*t^2)dt;
	std::complex<double> CalcM0_ERepulsion(const QuartetPara & qPara,const int m);
	std::complex<double> CalcM0_Overlap(const QuartetPara & qPara, const int m);
	std::complex<double> CalcM0_NA(const QuartetPara & qPara, const int m);
	void 				 CalcQuartetPara(const std::vector<PrimitiveGauss * > & ABCD, QuartetPara & qparameter ,const std::string & type);
	std::complex<double> CalcRM(const QuartetPara & qPara,const RM & rm, const std::string & type); //calculate the value [r]^m;
	std::complex<double> RMKernel(const QuartetPara & qPara,const RM & rm, const std::string & type); //calculate the value [r]^m;
	void 				 CalcABInt(std::vector<PrimitiveGauss *> CD, std::vector<std::complex<double> >::iterator & beginIter, std::vector<std::complex<double> >::iterator & endIter, const size_t stride);
	void 				 MDTransform1v(const MDTransformInfo &);
	void 				 MDTransform2v(const MDTransformInfo & mdTransInfo);
	void				 MDTransformVector(const std::vector<PrimitiveGauss *> &, std::vector<std::complex<double> >::iterator & beginIter, std::vector<std::complex<double> >::iterator & endIter, const size_t stride);
	void				 MDTransformTable(const std::vector<PrimitiveGauss *> &, std::vector<std::complex<double> > & pqTable);
	void 				 GetFullPQTable(const std::vector<PrimitiveGauss *> &ABCD, std::vector<std::complex<double> > & pqTable,const std::string & type);
	void 				 GetParIdx(const int & spaceBasis_i,int &parIdx_i,int & rel_idx);
	std::complex<double> GetMatrixElementFromTable(const std::vector<int> &, const std::vector<int> &, const std::vector<std::complex<double> > &);
	//*********for contract gauss basis**************//
	void 				 ParseGaussianBasis();
	void 				 ParseMomentum(const std::string & str, std::vector<int> & L);
	void 				 AddNewBasis(const std::vector<int> & L,const std::vector<std::complex<double> > & alpha, const std::vector<std::vector<double> > & coeff, double R[3]);
	void				 AdjustBasisSize(const std::vector<int> & L);
	void 				 CalcPrimitiveGaussNorm();
protected:
	int  CalcERepulsion(int piecen, std::vector<std::complex<double> > & rslt);
	std::complex<double> CalcABCDInt(const std::vector<PrimitiveGauss * > & ABCD, const std::string & type);//may be this is enough
};

#endif
