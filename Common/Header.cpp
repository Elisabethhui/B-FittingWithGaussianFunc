#ifndef _HEADER_CPP_
#define _HEADER_CPP_
#include "Header.h"

std::string trim(const std::string & str){//trim string, erase the tail and head space, and also \r\n at the end
	std::string whitespaces (" \t\f\v\n\r ");
	size_t n1=str.find_first_not_of(whitespaces);
	if(n1==str.npos){
		return "";//empty string
	}else{
		return str.substr(n1,str.find_last_not_of(whitespaces)-n1+1);
	}
}
#endif

