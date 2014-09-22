#ifndef _IO_H_
#define	_IO_H_
/*
IO is to control input and output flow.
 */

#include"../Parser/InputParser.h"
#include "./OutputController.h"
#include"../Parallelization/Parallelization.hpp"

class IO{
public:
	InputParser inpParser;//parser to analyze variable in input file.
	OutputController outPut;
	std::string inputFile;
	Parallelization * parallel;
//OutputM
private:
public://method
	int Init(std::string & inputfile,Parallelization * parallel);
	int PPrintf(const char * format, ...);//Parallel printf, only output through the process rank=0
	void Finalize();
private://method
};

#endif

