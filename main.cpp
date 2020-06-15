//	JUST FOR TESTING/DEBUGGING

#include <iostream>
using namespace std;

#include "includes/matrix.h"

typedef matrix<int> MI;

int main() {
	
	matrix<int> m=reshape(linspace(1,3),matrix<int>(2,1,{1,3}));
	
	try{
		cout<<m<<endl<<fft(m);
	}
	catch(const string& e){
		cout<<e;
	}

	return 0;
}

