#include"IO.hpp"

int IO::Init(std::string & inpFile,Parallelization * parallel){
	this->parallel=parallel;
	inputFile=inpFile;
	if(inpParser.Init(inpFile)) return 1;
	if(outPut.Init()) return 0;
	return 0;
}

void IO::Finalize(){
		outPut.Finalize();	
}

int IO::PPrintf(const char * format, ...){
	if(parallel->rank==0){
		va_list args;
		va_start(args, format);
		vprintf(format,args);
		va_end(args);
	}
	return 0;
}
