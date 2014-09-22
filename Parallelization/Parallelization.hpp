#ifndef PARALLELIZATION_HPP_
#define PARALLELIZATION_HPP_
#include<mpi.h>
#include<vector>
#include<cstring>
#include<complex>

class Parallelization{
public:
	MPI_Comm comm_world;
	bool Init(int argc, char ** argv);//Pass necessary arguments for MPI initialization
	bool MPIInit();
	bool Finalize();
	int AssembleData(void * buffer2Send,void * AssembleBuffer,int unitSize,const std::vector<int> & indexRange);//each process own a piece of data; this function make each process own a copy of whole data
	int AddCplxVecFromEachProc(std::complex<double> * vec2send, std::complex<double> * vec2receive,int size, int root);//Sum the vec from each process at processor with rank
	~Parallelization();
private:
public:
	int argc; // used for MPI initialization
	char ** argv;// used for MPI initialization
	int rank; //rank id of current processor
	int totalProc; //total number of rank.
};
#endif /* PARALLELIZATION_HPP_ */
