/*
The class InputParser is used to open a input file which contains all value for input variables.
The input variable can have three four kind of format like:
1.
variable=value;
 	value can be type of string, int, float, complex or mathematical expression.
2.

%variable
cl | c2 | c3 | c4
c1 | c2 | c3 | c4
%
block variable, column value c1, c2, c3, c4 can be any type mentioned above.
Anything after # will be omitted.

The following shows a example to get the valve of variable of each type.
Example1.
inp file contain the following variable
	length=10;
	boxsize=length;
To obtain the variable boxsize's value, just do
	obj.GetVariableValue("boxsize","double",&boxsize);

Example2.
%species
'H' | 1 | 'spec_ps_psf' | 1 | 0 | 0 |0
%
block  variable's value is treated as a table. To obtain a table item, just use
obj.GetBlockValue("blockVariableName",row,col,"typeName",&rslt);

Notice:
before query any variable, all classes should register it own reserved variable name
to avoid any mixing control variable.
bool IfReserved(variableName); check if variable is reserved;
ReserveVar(variableName, className); reserve it own variable, each variable can only reserved once.
Before getting variable value, use
IfVariableExist(variableName) to check if variable provided in input file.
For complex value cannot be nested. it is better to use block replace complex value.
Init() should be called before any other member function is called.
 */
#ifndef _INPUTPARSER_H_
#define	_INPUTPARSER_H_
#include<cstdio>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<list>
#include<complex>
#include<map>
#include<cstdlib>
#include<algorithm>
#include"./MathExprParser/parser.h"

struct ParItem{
	std::string name; //name of the paramters
	std::string type; //type value of the parameters;
	std::vector < std::vector < std::string > > value; //it is a string table
	void Clear();
};
class InputParser{
private:
	std::string filename;
	std::vector<ParItem> parameters;
	std::map<std::string, int> idxMap;
private:
	Parser prs;
	std::map<std::string, std::string> reservedVar;
	//store the variable name reserved by each class.
public: //method
	//construct
	InputParser(){};
	bool ReservedBy(const std::string & bar, std::string & name);
	int	 Init(const std::string fileName);
	bool IfVariableExist(const std::string &) const;
	bool IfReserved(const std::string &);
	void ReserveVar(const std::string & var,const std::string & className);
	void Print(std::ofstream * );
	void Print(const std::string &);
	void GetVariableValue(const std::string & val, const std::string & type, void * rslt);
	void GetBlockValue(const std::string & var, const int x,const int y, const std::string & type, void * rslt);
	void GetBlockSize(const std::string & var, int & x, int & y);
	void GetBlockRowLen(const std::string & var, int & row);
	void GetBlockColLen(const std::string & var, const int row, int & col);
	int  IsBlockVariable(const std::string & var);
	bool IsStringVariable(const std::string & str);
	//evaluate the string expression and return the their value.
private: //method
	void Add2IdxMap(const std::string &);
	void ParserInit();
	int String2TypeX(const std::string & expr, const std::string & type, void * rslt);
};

#endif
