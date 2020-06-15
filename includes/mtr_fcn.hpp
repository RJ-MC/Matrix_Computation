#ifndef mtr_fcn_hpp
#define mtr_fcn_hpp

#include "matrix.h"

#include <cmath>

template <class _T>
matrix<_T> chs(const matrix<_T>& m, bool (*keep_fcn)(const _T&,const _T&)){
	auto r=m.get(linspace(0,(size(m).get(1))-1)*size(m).get(0));
	for(int i=1;i<size(m).get(0);++i)
		for(int j=0;j<size(m).get(1);++j)
			r[0][j]=keep_fcn(r[0][j],m.get(i,j))?r[0][j]:m.get(i,j);
	return r;
}
template <class _T>
matrix<_T> max(const matrix<_T>& m){
	return chs(m,(bool(*)(const _T&,const _T&))([](const _T& a,const _T& b)->bool{return a>=b;}));
}
template <class _T>
matrix<_T> min(const matrix<_T>& m){
	return chs(m,(bool(*)(const _T&,const _T&))([](const _T& a,const _T& b)->bool{return a<=b;}));
}


template <class _T>
auto vecnorm(const matrix<_T>& m, double p=2) -> matrix<typename remove_const<decltype(pow(_T(1),1.5))>::type>& {
	auto n=new matrix<typename remove_const<decltype(pow(_T(1),1.5))>::type>(1,size(m).get(1));
	auto fst_row=n->operator[](0);
	for (int i = 0; i < size(m).get(0); ++i)
		for (int j = 0; j < size(m).get(1); ++j)
			fst_row[j] += pow(abs(m.get(i,j)),p);
	for (int j = 0; j < size(m).get(1); ++j)
		fst_row[j]=pow(fst_row[j],1.0/p);
	return *n;
}

template <class _T>
matrix<_T>& diag(const matrix<_T>& m,int n=0){
	auto d=new matrix<_T>(1,MAX(0,(n<0)?(MIN(size(m).get(0),size(m).get(1)+n)):(MIN(size(m).get(1),size(m).get(0)-n))));
	if(n>=0)
		for(int i=0;i<size(*d).get(1);++i)
			(*d)[0][i]=m.get(i+n,i);
	else
		for(int i=0;i<size(*d).get(1);++i)
			(*d)[0][i]=m.get(i,i-n);
	return *d;
}

template <class _T>
matrix<_T> utri(const matrix<_T>& m,int n=0){
	auto u=new matrix<_T>(m);
	int r(MAX(0,n)),c(MAX(0,-n)),curr_c;
	while(r<size(m).get(0)){
		curr_c=0;
		while(curr_c<r-n)
			(*u)[r][curr_c++]=0;
		++r;
	}
	return *u;
}

template <class _T>
matrix<_T> ltri(const matrix<_T>& m,int n=0){
	auto l=new matrix<_T>(m);
	int r(0),curr_c;
	while(r<size(m).get(0)&&(r-n+1)<size(m).get(1)){
		curr_c=r-n+1;
		while(curr_c<size(m).get(1))
			(*l)[r][curr_c++]=0;
		++r;
	}
	return *l;
}

#endif
