#include <iostream>
using namespace std;

#include "includes/matrix.h"

typedef matrix<int> MI;

int main() {
	
	auto m=reshape(linspace(1,12),MI(2,1,{3,4}));
	
	try{
		cout<<m<<ltri(m,1);
	}
	catch(const string& e){
		cout<<e;
	}

	return 0;
}

