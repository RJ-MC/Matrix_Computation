#ifndef mtr_fcn_hpp
#define mtr_fcn_hpp

#include <cmath>
#include <complex>

#include "CONSTANTS.h"
#include "matrix.h"

// Trick to allow type promotion
template <typename T>
struct identity_t { typedef T type; };

// Make std::complex<> work better....
#define COMPLEX_OPS(OP)								\
template <typename _Tp>								\
std::complex<_Tp>									\
operator OP(std::complex<_Tp> lhs, 					\
		const typename identity_t<_Tp>::type & rhs){\
    return lhs OP rhs;								\
}													\
template <typename _Tp>								\
std::complex<_Tp>									\
operator OP(const typename identity_t<_Tp>::type 	\
			& lhs, const std::complex<_Tp> & rhs){	\
	return lhs OP rhs;								\
}
COMPLEX_OPS(+)
COMPLEX_OPS(-)
COMPLEX_OPS(*)
COMPLEX_OPS(/)
#undef COMPLEX_OPS

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
matrix<_T>& utri(const matrix<_T>& m,int n=0){
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
matrix<_T>& ltri(const matrix<_T>& m,int n=0){
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

template <class _T>
matrix<_T>& sum(const matrix<_T>& m){
	matrix<_T>* r;
	if(size(m).get(0)==1){
		r=new matrix<_T>(1,1);
		for(int i=0;i<size(m).get(1);++i)
			(*r)[0][0]+=m.get(0,i);
	}else{
		r=new matrix<_T>(m.get(linspace(0,size(m).get(1)-1)*size(m).get(0)));
		for(int i=1;i<size(m).get(0);++i)
			for(int j=0;j<size(m).get(1);++j)
				(*r)[0][j]+=m.get(i,j);
	}
	return *r;
}

template <class _T>
matrix<_T>& cmb(const matrix<_T>& m1, const matrix<_T>& m2, bool _horizontal=true){
	matrix<_T>* r;
	if(_horizontal){
		if(size(m1).get(0)!=size(m2).get(0))
			throw(string("[Error] cmb : Number of rows should be same."));
		r=new matrix<_T>(size(m1).get(0),size(m1).get(1)+size(m2).get(1));
		r->set(linspace(0,numel(m1)-1),m1);
		r->set(linspace(0,numel(m2)-1)+numel(m1),m2);
	}else{
		if(size(m1).get(1)!=size(m2).get(1))
			throw(string("[Error] cmb : Number of columns should be same."));
		r=new matrix<_T>(size(m1).get(0)+size(m2).get(0),size(m1).get(1));
		r->set(linspace(0,size(m1).get(0)-1),linspace(0,size(m1).get(1)-1),m1);
		r->set(linspace(0,size(m2).get(0)-1)+size(m1).get(0),linspace(0,size(m1).get(1)-1),m1);
	}
	return *r;
}

template <class _T>
auto fft2(const matrix<_T>& m)->matrix<typename remove_const<decltype(complex<double>(1)*_T(1))>::type>&{
	int row(size(m).get(0)),col(size(m).get(1));
	if(row/2*2<row)
		throw(string("[Error] fft2 : Number of rows should be a power of 2."));
	auto r=new matrix<typename remove_const<decltype(complex<double>(1)*_T(1))>::type>(m);
	complex<double> w=exp(complex<double>(0,-1)*PI*2.0/row);
	if(row>1){
		auto yt=fft2(m.get(linspace(0,row-1,2),linspace(0,col-1)));
		auto yb=fft2(m.get(linspace(1,row-1,2),linspace(0,col-1)));
		auto z=((w^linspace(0,row/2-1)).T(0)^ones(matrix<int>(2,1,{1,col})))*yb;
		r->set(linspace(0,row/2-1),linspace(0,col-1),yt+z);
		r->set(linspace(row/2,row-1),linspace(0,col-1),yt-z);
	}
	return *r;
}



#endif
