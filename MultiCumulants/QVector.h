#ifndef MULTICUMULANTS_QVECTOR_H
#define MULTICUMULANTS_QVECTOR_H

#include <string>

class QVector
{
public:
	virtual const char* name() const { return "Q Vector"; }
	virtual const char* classname() const {  return "QVector"; }
	QVector() {}
	~QVector() {}


	std::string toString(){
		std::string s = "";
		s += classname();
		s += "<>";
		return s;
	}
	
};



#endif