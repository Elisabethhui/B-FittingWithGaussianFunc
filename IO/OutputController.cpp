#ifndef _OUTPUTCONTROLLER_CC_
#define _OUTPUTCONTROLLER_CC_

#include "./OutputController.h"

bool OutputController::IfOutputFileExist(const std::string & file)
{
	std::map<std::string, std::ofstream * >::iterator iter=outFileList.find(file);
	return	(iter!=outFileList.end())? true:false;
}



bool OutputController::CloseFile(const std::string & file)
{
	std::map<std::string,std::ofstream *>::iterator iter=outFileList.find(file);
	if(iter!=outFileList.end()){//file exist
		(*iter).second->close();
		delete (*iter).second;
		outFileList.erase(iter);
		return true;
	}else{
		printf("warning::trying to close file %s which does not exist\n",file.c_str());
	return true;
	}
}

std::ofstream * OutputController::GetOutputFile(const std::string & file)
{//if file exist, return the ofstream, otherwise NULL:
	std::map<std::string, std::ofstream * >::iterator iter=outFileList.find(file);
	return (iter!=outFileList.end())? (*iter).second:NULL;
}

std::ofstream	* openFile(const std::string & file)
{	std::ofstream * ofile=NULL;
	std::string filePath=trim(file);
	if(filePath.substr(0,2).compare("./")==0){//remove "./";
			filePath=filePath.substr(2,filePath.size()-2);
	}
	size_t lastPos=filePath.find_last_of('/');
	if(lastPos!=filePath.npos){
		boost::filesystem::path dir(filePath.substr(0,lastPos));
		if(!boost::filesystem::exists(dir)){
			if(!boost::filesystem::create_directories(dir)){
				std::cerr<<"fail to create folder"<<std::endl;
			}
		}
	}

	ofile=new std::ofstream(); //need to be delete in destructor.
	(*ofile).open(trim(file).c_str(),std::ofstream::out);
	if(!(*ofile).good()){
		std::cerr<<"fail to create file "<<file<<std::endl;
		exit(1);
	}
	return ofile;
}
/*
std::ofstream * openFile2(const std::string & file)
{//open a file with w mode. Support open file in subdirectory.
	std::ofstream * ofile=NULL;
	std::string filePath=trim(file);
	if(filePath.substr(0,2).compare("./")==0){
		filePath=filePath.substr(2,filePath.size()-2);
	}
	size_t pos;
	size_t startPos=0;
	std::string dir;
	size_t lastPos=filePath.find_last_of('/');
	pos=filePath.find_first_of('/');
	if(lastPos!=filePath.npos){
		while(pos<=lastPos){
			dir=filePath.substr(0,pos);
			if(mkdir(dir.c_str())!=0){
				std::cerr<<"Fail to create directory "<<dir<<": ";
				std::cerr<<strerror(errno)<<std::endl;
			//	exit(1);
			}
			startPos=pos+1;
			pos=filePath.find_first_of('/',startPos);
		}
	}
	ofile=new std::ofstream();
	(*ofile).open(trim(file).c_str(),std::ofstream::out);
	if(!(*ofile).good()){
		std::cerr<<"fail to create file "<<file<<std::endl;
		exit(1);
	}
	return ofile;
}
*/

std::ofstream * OutputController::AddOutputFile(const std::string & file){
	std::ofstream * ofile=GetOutputFile(file);
	if(NULL!=ofile){
		std::cout<<"Warning: output file already exist"<<std::endl;
		return ofile;
	}else{
		std::pair<std::string, std::ofstream * > outFilePair(file,openFile(file));
		outFileList.insert(outFilePair);
		return outFilePair.second;
	}
}

std::ofstream * OutputController::AddOutputFile(const std::string & name,const std::string &path){
	std::ofstream * ofile=GetOutputFile(name);
	if(NULL!=ofile){
		printf("Warning: output file %s already exist\n",path.c_str());
		return ofile;
	}else{
		std::pair<std::string, std::ofstream *> outFilePair(name,openFile(path));
		outFileList.insert(outFilePair);
		return outFilePair.second;
	}
}

int OutputController::Init()
{//open some file for certain type output.
/*	std::string parserInfo="./Parser/parser";
	std::string eigen="./Eigen/EigenValue";
	std::string eigenVec="./Eigen/EigenVec";
	AddOutputFile(parserInfo);
	AddOutputFile(eigen);
	AddOutputFile(eigenVec);*/
	return 0;
}

OutputController::~OutputController()
{//close all opened file

}


void OutputController::Finalize(){
	std::map<std::string, std::ofstream *>::iterator iter;
	for(iter=outFileList.begin();iter!=outFileList.end();iter++)
	{
		try{
			(*iter).second->close();
			delete (*iter).second;
		}catch(std::exception  & e ){
			std::cerr<<"fail to close file "<<(*iter).first<<" "<<e.what()<<std::endl;
		}
	}

}

#endif

