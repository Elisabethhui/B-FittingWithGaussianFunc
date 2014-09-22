#ifndef _OUTPUT_CONTROLLER_H_
#define _OUTPUT_CONTROLLER_H_
#include<iostream>
#include<fstream>
#include<string>
#include<map>
#include<vector>
#include<cstring>
#include"../Common/Header.h"
#include<cstdlib>
#include<boost/filesystem.hpp>

class OutputController{
public://methods
		~OutputController();//all output file will be closed
		std::ofstream * AddOutputFile(const std::string &);
		std::ofstream * GetOutputFile( const std::string & );
		std::ofstream *	AddOutputFile(const std::string & name,const std::string &path);
		bool CloseFile(const std::string & path);
		bool IfOutputFileExist(const std::string &);
		int	 Init();
		void Finalize();
private://methods
public://properties
private://properties
	std::map<std::string, std::ofstream * > outFileList;
};

#endif
