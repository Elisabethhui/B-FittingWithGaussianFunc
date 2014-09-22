#ifndef _INPUTPARSER_
#define	_INPUTPARSER_
#include"InputParser.h"
#include "../Common/Header.h"
#include<boost/filesystem.hpp>

typedef std::pair<std::string, int> IdxPair;
void InputParser::Add2IdxMap(const std::string & varName)
{//add variable index to hash table
// if variable exist, then emitt error and exit
	if(idxMap.find(varName)!=idxMap.end()){
		std::cerr<<"redefined: "<<varName<<std::endl;
		exit(1);
	}
	IdxPair idxpair(varName,parameters.size()-1);
	idxMap.insert(idxpair);
}

void ParItem::Clear(){
	this->name="";
	this->type="";
	this->value.clear();
}
int InputParser::Init(const std::string fileName)
{//function init() is to open file and read all input variable to the list
//parameters. All value is stored as string.
	filename=fileName;
	if(!boost::filesystem::exists(filename)){
		std::cout<<filename<<" does not exist, please check your input file"<<std::endl;
		return 1;
	}

	std::ifstream file(filename.c_str(),std::ifstream::in);
	std::string line;
	ParItem parItem;
	bool ifInBlock=false;
	size_t n;
	while(file.good()){
		std::getline(file,line); //delimiter here is to omit everything after #
		line=trim(line); //remove extra space and \r delimeter
		if((n=line.find('#'))!=line.npos) line=line.substr(0,n);
		if(line.empty()) continue;//if empty, then ,just go to next line.
		if('%'==line[0]){
			if(ifInBlock)
			{
				if(1!=line.size()){
				//IO??
					std::cout<<line<<std::endl;
					perror(("Please check your block format: "+parItem.name).c_str());
					return 1; // MPI terminate signal has to be emitted.
				}else{//end of block
					parameters.push_back(parItem);
					Add2IdxMap(parItem.name);
					parItem.Clear();
					ifInBlock=false;
					continue;
				}
			}else{
				if(1==line.size()){
					//IO??
					printf("Block name has to be specified\n");
					return 1;
				}else{
					parItem.name=trim(line.substr(1,line.size()-1));
					parItem.type="block";
					ifInBlock=true;
					continue;
				}
			}
		}else{
			if(ifInBlock){//a new block line
				std::istringstream istr(line);
				std::string col;
				std::vector<std::string> blockLine;
				while(istr.good())
				{
					getline(istr,col,'|');//get column value
					if(col.empty()){
					//IO??
						perror(("Please check your column value of block: "+parItem.name).c_str());
						return 1;
					}else{
						blockLine.push_back(trim(col));
					}
				}
				parItem.value.push_back(blockLine);
				blockLine.clear();
			}else{// a new non-block variable;
				size_t n=line.find('=');
				if(std::string::npos==n || n==line.size()-1){
					//IO??
					printf("variable should be assigned\n");
					return 1;
				}
				if(0==n){
					//IO??
					printf("variable name is missing\n");
					return 1;
				}
				parItem.name=trim(line.substr(0,n));
				parItem.type="variable";
				std::vector<std::string> vec;
				vec.push_back(trim(line.substr(n+1,line.size()-n-1)));
				parItem.value.push_back(vec);
				parameters.push_back(parItem);
				Add2IdxMap(parItem.name);
				parItem.Clear();
			}
		}
	}
	file.close();
	ParserInit();
	return 0;
}

void InputParser::Print(std::ofstream * ofile)
{
	if((*ofile).good()){
		for(std::vector<ParItem>::iterator iter=parameters.begin();iter!=parameters.end();iter++)
		{//output variable name and their value
			if((*iter).type!="block")
			{//nonblock variable
				(*ofile)<<(*iter).name<<"="<<(*iter).value.front().front()<<std::endl;
			}else
			{//block variable
				(*ofile)<<(*iter).name<<"=";
				for(size_t n=0; n<(*iter).value.size(); n++){
					(*ofile)<<std::endl<<"\t";
					for(size_t m=0; m<(*iter).value[n].size();m++){
						(*ofile)<<(*iter).value[n][m]<<"|";
					}
				}
				(*ofile)<<std::endl;
			}
		}
	}else{
		std::cerr<<"input parameter: ofstream ofile is not valid"<<std::endl;
	}
}

void InputParser::Print(const std::string & filePath)
{//print all parameters information.
	std::ofstream ofile(filePath.c_str(),std::ofstream::out);
	if(ofile.good()){
		Print(&ofile);
		ofile.close();
	}else{
		std::cerr<<"InputParser::Fail to open file "<<filePath<<std::endl;
	}

}

bool InputParser::IsStringVariable(const std::string & str)
{//if string contain single quote or double quote, then it should be string variable
//or if the variable name does not defined, it is also treated as a string
	if(str.find("\"") !=std::string::npos || str.find("'") !=std::string::npos) return true;
	return false;
}
bool IsComplex(const std::string & str,std::string & real,std::string & imaginary)
{//if it is complex, return true and real part and imaginary part
//if it is not complex, return false
//complex format (real, imaginary)
	std::string t_str=trim(str);
	size_t pos;
	if(t_str[0]=='(' && t_str[t_str.size()-1]==')' && (pos=t_str.find(','))!=t_str.npos){
		real=t_str.substr(1,pos-1);
		imaginary=t_str.substr(pos+1,t_str.size()-pos-2);
		return true;
	}else{
		return false;
	}
}

int InputParser::IsBlockVariable(const std::string & var){
//return value 1 is block variable
//if not exist or not block variable, 0 is returned
	std::map<std::string, int>::iterator iter;
	iter=idxMap.find(var);
	if(iter!=idxMap.end()){
		if(parameters[(*iter).second].type.compare("block")==0){
			return 1;
		}
	}
	return 0;
}

void InputParser::ParserInit()
{//for non-block variable and non-string variable, evaluate them first.
	char * rslt=NULL;
	std::string *pStr=NULL;
	std::string real,imaginary;
	int flag=0;
	for(std::vector< ParItem >::iterator iter=parameters.begin(); iter!=parameters.end();iter++)
	{
		if((*iter).type=="block"){
			for(size_t i=0;i<(*iter).value.size();i++)
				for(size_t j=0;j<(*iter).value[i].size();j++)
				{//check each block column, if they are math expression
					pStr=&((*iter).value[i][j]);
					if(IsStringVariable(*pStr)) continue;
					if(IsComplex(*pStr,real,imaginary))
					{//complex
						real=prs.parse(real.c_str());
						imaginary=prs.parse(imaginary.c_str());
						(*pStr)="("+real+","+imaginary+")";
					}else{//not complex
						rslt=prs.parse((*pStr).c_str());
						(*pStr)=rslt;
					}
				}
		}else{
			pStr=&((*iter).value[0][0]);
			if(IsStringVariable(*pStr)) continue; //if it is string variable, also pass
			if(IsComplex(*pStr,real,imaginary))
			{//complex
				real=prs.parse(real.c_str());
				imaginary=prs.parse(imaginary.c_str());
				(*pStr)="("+real+","+imaginary+")";
			}else
			{//not complex
				rslt=prs.parse(((*iter).name+"="+*pStr).c_str());
				*pStr=rslt;
			}
		}
	}
}

bool InputParser::ReservedBy(const std::string & var, std::string & name){
	std::map<std::string, std::string>::iterator iter;
	iter=reservedVar.find(var);
	if(iter!=reservedVar.end()){
		name=(*iter).second;
		return true;
	}else{
		return false;
	}
}
bool InputParser::IfReserved(const std::string & var)
{//check if the variable name is reserved by other classes.
		std::map<std::string, std::string>::iterator iter;
		iter=reservedVar.find(var);
		if(iter!=reservedVar.end()){
			std::cout<<var<<" reserved by class: "<<(*iter).second<<std::endl;
			return true;
		}
		return false;
}

void InputParser::ReserveVar(const std::string & var, const std::string & className)
{//put a reserved variable in reservedVar;
	std::string name;
	if(ReservedBy(var,name))
	{//var is reserved
		if(className.compare(name)!=0){
			throw "variable is already reserved";
		}
	}else
	{//not reserved, then register it
		std::pair<std::string, std::string> resItemPair(var,className);
		reservedVar.insert(resItemPair);
	}
}
bool InputParser::IfVariableExist(const std::string & var) const

{//check if variable exist in input file
	if(idxMap.find(var)!=idxMap.end()){
		return true;
	}else{
		return false;
	}
}
int InputParser::String2TypeX(const std::string & expr, const std::string & type, void * rslt)
{
	if("int"==type){
		int * iRslt=static_cast<int *>(rslt);
		 (*iRslt)=strtol(expr.c_str(),NULL,10);
	}else if("double"==type){
		double * fRslt=static_cast<double *>(rslt);
		(*fRslt)=strtod(expr.c_str(),NULL);
	}else if("string"==type){
		std::string * strRslt=static_cast<std::string *>(rslt);
		size_t x1,x2;
		x1=expr.find_first_of("\"'");
		x2=expr.find_last_of("\"'");
		if(x1==x2){
			std::cout<<"string format is not correct: "<<expr<<std::endl;
			exit(1);
		}
		(*strRslt)=expr.substr(x1+1,x2-x1-1);//need to eliminate the first and last quote
	}else if("complex"==type){
		std::complex<double> * cplxRslt=static_cast<std::complex< double > * >(rslt);
		size_t pos=expr.find(",");
		(*cplxRslt)=std::complex<double>(strtod(expr.substr(1,pos-1).c_str(),NULL),
							strtod(expr.substr(pos+1,expr.size()-pos-2).c_str(),NULL));
	}else if("bool"==type){
		bool * bRslt=static_cast<bool *>(rslt);
		size_t x1,x2;
		x1=expr.find_first_not_of("\"");
		x2=expr.find_last_not_of("\"");
		std::string t_expr(expr.substr(x1,x2-x1+1));
		std::transform(t_expr.begin(),t_expr.end(),t_expr.begin(),::tolower);
		(*bRslt)=(0==t_expr.compare("yes")||0==t_expr.compare("true"))? true:false;
	}else{
		std::cout<<"unsupported type: "<<type<<std::endl;
		return 1;
	}
	return 0;
}
void InputParser::GetVariableValue(const std::string & var, const std::string & type, void * rslt)
{//get the value of nonBlock variable.
//if no error happen, return 0;
	std::string & expr=(parameters[idxMap[var]].value[0][0]);
	if(1==String2TypeX(expr,type,rslt))
		exit(1);
}
void InputParser::GetBlockValue(const std::string &var, const int x, const int y, const::string &type,void * rslt)
{//get block value at row x and coloumn y.
//if no error, return 0;
	try{
		std::string & expr=parameters[idxMap[var]].value[x-1][y-1];
		if(1==String2TypeX(expr,type,rslt))
			exit(1);
	}
	catch(std::exception & e){
		std::cout<<e.what()<<std::endl;
		std::cout<<"unexpected error when trying to get block value at row "<<x<<" col "<<y<<std::endl;
	}
}

void InputParser::GetBlockSize(const string & var, int & x, int & y)
{
	ParItem * pParItem=&parameters[idxMap[var]];
	x=pParItem->value.size();
	y=pParItem->value.front().size();
}

void InputParser::GetBlockRowLen(const string & var, int & row){
	ParItem * pParItem=&parameters[idxMap[var]];
	row=pParItem->value.size();
}

void InputParser::GetBlockColLen(const string & var, const int row, int & col){
//get the size at row
		ParItem & parItem(parameters[idxMap[var]]);
		col=parItem.value[row-1].size();//the index for block is start from 1;
}


#endif

