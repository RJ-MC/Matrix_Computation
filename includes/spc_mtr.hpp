#ifndef spc_mtr_hpp
#define spc_mtr_hpp

#include "matrix.h"

matrix<int> eye(int n) {
	matrix<int> I(n, n);
	for (int i = 0; i < n; ++i)
		I.mtr[i][i] = 1;
	return I;
}

matrix<int> zeros(const matrix<int> &size) {
	matrix<int> I(size,{});
	return I;
}

matrix<int> ones(const matrix<int> &size) {
	matrix<int> I(size,{});
	for (int i = 0; i < size.get(0); ++i)
		for(int j=0;j<size.get(1);++j)
			I.mtr[i][j] = 1;
	return I;
}

template <class _U,class _V,class _W=int>
auto linspace(const _U& start,const _V& end,const _W& step=1)->matrix<typename remove_const<decltype(_U(0)+_W(0))>::type>{
	matrix<typename remove_const<decltype(_U(0)+_W(0))>::type> r(1,(int)((end-start)/step)+1);
	decltype(_U(0)+_W(0)) curr(start);
	for(int i=0;i<size(r).get(1);++i,curr+=step)
		(r)[0][i]=curr;
	return r;
}

template <class _U>
matrix<_U> linspace(const _U& start){
	return matrix<_U>(1,1,{start});
}

#endif
