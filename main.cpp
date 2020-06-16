//	JUST FOR TESTING/DEBUGGING

#include <iostream>
using namespace std;

#include "includes/matrix.h"

typedef matrix<int> MI;

template <class _T>
matrix<_T>& fht(const matrix<_T>& m){
	int row(size(m).get(0)),col(size(m).get(1));
	if(row/2*2<row)
		throw(string("[Error] fft2 : Number of rows should be a power of 2."));
	auto r=new matrix<_T>(m);
	auto x=m.get(linspace(row/2,row-1),linspace(0,col-1));
	if(row>1){
		auto z=fht(m.get(linspace(0,row/2-1),linspace(0,col-1)));
		r->set(linspace(0,row-1,2),linspace(0,col-1),z+x);
		r->set(linspace(1,row-1,2),linspace(0,col-1),z-x);
	}
	return *r;
}

int main() {
	
	matrix<int> m=reshape(linspace(1,24),matrix<int>(2,1,{8,3}));
	
	try{
		cout<<m<<endl<<fht(m);
	}
	catch(const string& e){
		cout<<e;
	}

	return 0;
}

