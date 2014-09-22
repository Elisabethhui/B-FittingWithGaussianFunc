#include"./Parallelization.hpp"
#include<complex>

bool Parallelization::Init(int argc, char ** argv)//Pass necessary arguments for MPI initialization
{
	this->argc=argc;
	this->argv=argv;
	return MPIInit();
}

bool Parallelization::MPIInit(){
	//init mpi
	int flag;
	MPI_Initialized(&flag);
	if( !flag && MPI_Init(&(this->argc),&(this->argv))!=0){
		printf("Parallelization::Init() fail to init MPI environment\n");
		return false;
	}
	MPI_Comm_rank(MPI_COMM_WORLD,&(this->rank));
	MPI_Comm_size(MPI_COMM_WORLD,&(this->totalProc));
}

bool Parallelization::Finalize(){
	int flag;
	MPI_Finalized(&flag);
	if(!flag){
		if(0!=MPI_Finalize()){
			printf("Parallelization::Finalize() fail to Finalize MPI environment\n");
			return false;
		}
	}
	return true;
}

int Parallelization::AssembleData(void * buffer2Send,void * assembleBuffer,int unitSize,const std::vector<int> & indexRange)
{//this function assemble the data stored in buffer2Send at each process to assembleBuffer
//indexRange is the global data index for each buffer2send in assembleBuffer
	void * p=assembleBuffer;
	MPI_Barrier(MPI_COMM_WORLD);
	for(size_t i=0;i<totalProc;++i){
		if(i==rank){
			MPI_Bcast(buffer2Send,(indexRange[i+1]-indexRange[i])*unitSize,MPI_BYTE,i,MPI_COMM_WORLD);
			memcpy(p,buffer2Send,unitSize*(indexRange[i+1]-indexRange[i]));
		}else{
			MPI_Bcast(p,(indexRange[i+1]-indexRange[i])*unitSize,MPI_BYTE,i,MPI_COMM_WORLD);
		}
		p+=unitSize*(indexRange[i+1]-indexRange[i]);
	}
	return 0;
}

Parallelization::~Parallelization()
{
	Finalize();
}

int Parallelization::AddCplxVecFromEachProc(std::complex<double> * vec2send, std::complex<double> * vec2receive, int size, int root){
//sum the complex vec from each proc and store it at processor with rank
	MPI_Barrier(MPI_COMM_WORLD);
	int sendcount=size*sizeof(std::complex<double>);
	std::complex<double> * pbuff(NULL);
	int receivecount=0;
	std::vector<std::complex<double> > buff;
	if(rank==root){
		buff.resize(size*totalProc);
		pbuff=&buff[0];
		receivecount=sendcount;
	}
	
	MPI_Gather( static_cast<void *>(vec2send),sendcount,MPI_BYTE,static_cast<void *>(pbuff),receivecount, MPI_BYTE,root,MPI_COMM_WORLD);

	if(rank==root){//add vec
		for(int n=0;n<size;++n){
			vec2receive[n]=0;
			for(int pn=0;pn<totalProc;++pn){
				vec2receive[n]+=pbuff[pn*size+n];
			}
		}
	}
	return 0;
}

