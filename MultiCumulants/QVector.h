#ifndef MULTICUMULANTS_QVECTOR_H
#define MULTICUMULANTS_QVECTOR_H

#include <string>

// logging library
#define LOGURU_IMPLEMENTATION 1
#include "loguru.hpp"

class QVector
{
public:
	virtual const char* name() const { return "Q Vector"; }
	virtual const char* classname() const {  return "QVector"; }
	QVector() {
		LOG_SCOPE_FUNCTION(INFO);
	}
	~QVector() {
		LOG_SCOPE_FUNCTION(INFO);
	}


	std::string toString(){
		std::string s = "";
		s += classname();
		s += "<>";
		return s;
	}
	
};



#endif