

#include "MultiCumulants/QVector.h"
#include "MultiCumulants/EtaRegion.h"

#include <iostream>
using namespace std;


int 
main(int argc, char** argv) {

	QVector qv;
	EtaRegion er(-1.5, -1.0);

	cout << qv.toString() << endl;
	cout << er.toString() << endl;

	return 0;
}