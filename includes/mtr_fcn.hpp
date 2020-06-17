#ifndef mtr_fcn_hpp
#define mtr_fcn_hpp

#include <cmath>
#include <complex>
#include <tuple>

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

template <class _T, class _U>
matrix<_U>& apply(const matrix<_T>& m, _U (*f)(const _T&)){
	matrix<_U>* r=new matrix<_U>(size(m),{});
	for(int i=1;i<size(m).get(0);++i)
		for(int j=0;j<size(m).get(1);++j)
			(*r)[i][j]=f(m.get(i,j));
	return *r;
}

matrix<int>& find(const matrix<bool>& m,bool _full_search=true){
	matrix<int> *r;
	if(sum(sum(matrix<int>(m))).get(0)){
		if(_full_search){
			r=new matrix<int>(1,sum(sum(matrix<int>(m))).get(0));
			int _curr=0;
			for(int i=0;i<numel(m);++i)
				if(m.get(i))
					r->set(linspace(_curr++),linspace(i));
		}else{
			int i=-1;
			while(!m.get(++i));
			r=new matrix<int>(1,1,{i});
		}	
	}else
		r=new matrix<int>(1,0);
	return *r;
}

template <class _T, class _V>
matrix<typename remove_const<decltype(_T(0)+_V(0))>::type>& cmb(const matrix<_T>& m1, const matrix<_V>& m2, bool _horizontal=true){
	matrix<typename remove_const<decltype(_T(0)+_V(0))>::type>* r;
	if(_horizontal){
		if(size(m1).get(0)!=size(m2).get(0))
			throw(string("[Error] cmb : Number of rows should be same."));
		r=new matrix<typename remove_const<decltype(_T(0)+_V(0))>::type>(size(m1).get(0),size(m1).get(1)+size(m2).get(1));
		r->set(linspace(0,numel(m1)-1),m1);
		r->set(linspace(0,numel(m2)-1)+numel(m1),m2);
	}else{
		if(size(m1).get(1)!=size(m2).get(1))
			throw(string("[Error] cmb : Number of columns should be same."));
		r=new matrix<typename remove_const<decltype(_T(0)+_V(0))>::type>(size(m1).get(0)+size(m2).get(0),size(m1).get(1));
		r->set(linspace(0,size(m1).get(0)-1),linspace(0,size(m1).get(1)-1),m1);
		r->set(linspace(0,size(m2).get(0)-1)+size(m1).get(0),linspace(0,size(m1).get(1)-1),m2);
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

template <class _T>
tuple<matrix<_T>, matrix<int>> simplex(const matrix<_T>& m, const matrix<int>& basic){
	auto M=new matrix<_T>(m);
	auto B=new matrix<int>(basic);
	//cout<<(*M);
	int row(size(m).get(0)),col(size(m).get(1));
	for(int i=0;i<numel(*B);++i)
		if(IS_NON_ZERO(M->get(row-1,B->get(i))))
			M->set(linspace(1,col)*row-1,
				M->get(linspace(1,col)*row-1)-M->get(linspace(0,col-1)*row+i)*M->get(row-1,B->get(i)));
	int _neg_col,_min_row;
	matrix<_T> _positive_row(0,0),_positive_row_val(0,0);
	while(sum(M->get(linspace(1,col-1)*row-1)<0).get(0)){
		_neg_col=find(M->get(linspace(1,col-1)*row-1)<0,false).get(0);
		if(!numel(find(M->get(linspace(0,row-2),linspace(_neg_col))>0,false)))
			throw(string("[Error] simplex : The problem is unbounded."));
		_positive_row=find(M->get(linspace(0,row-2),linspace(_neg_col))>0);
		_positive_row_val=M->get(_positive_row,linspace(col-1))/M->get(_positive_row,linspace(_neg_col));
		_min_row=_positive_row.get(find(_positive_row_val==min(_positive_row_val).get(0),false).get(0));
		M->set(linspace(0,col-1)*row+_min_row,
				M->get(linspace(0,col-1)*row+_min_row)/M->get(_min_row,_neg_col));
		for(int i=0;i<row;++i)
			if(i!=_min_row)
				M->set(linspace(0,col-1)*row+i,
					M->get(linspace(0,col-1)*row+i)-M->get(linspace(0,col-1)*row+_min_row)*M->get(i,_neg_col));
		B->set(matrix<int>(1,1,{_min_row}),matrix<int>(1,1,{_neg_col}));
	}
	return make_tuple(*M,*B);
}

template <class _T>
auto simplex(const matrix<_T>& m)->decltype(-1*_T(1)){
	auto M=new matrix<_T>(m);
	M->set(linspace(numel(*M)-1),linspace(0));
	int row(size(*M).get(0)),col(size(*M).get(1));
	auto _neg_row=find(M->get(linspace(0,row-2),linspace(col-1))<0);
	M->set(_neg_row,linspace(0,col-1),-1*M->get(_neg_row,linspace(0,col-1)));
	auto c=M->get(linspace(row-1),linspace(0,col-2));
	auto b=M->get(linspace(0,row-1),linspace(col-1));
	(*M)=cmb(cmb(cmb(M->get(linspace(0,row-2),linspace(0,col-2)),eye(row-1)),
			cmb(zeros(matrix<int>(2,1,{1,col-1})),ones(matrix<int>(2,1,{1,row-1}))),false),b);
	auto SOL=simplex((*M),linspace(col-1,col+row-3));
	(*M)=get<0>(SOL);
	auto B=get<1>(SOL);
	cout<<"End of phase I :\n"<<(*M);
	if(IS_NON_ZERO(M->get(numel(*M)-1)))
		throw(string("[Error] simplex : No feasible solution."));
	b=M->get(linspace(0,row-1),linspace(size(*M).get(1)-1));
	*M=cmb(cmb(M->get(linspace(0,row-2),linspace(0,col-2)),c,false),b);
	(*M)=get<0>(simplex((*M),B));
	cout<<"End of phase II :\n"<<(*M);
	return -1*(*M).get(numel(*M)-1);
}

template<class _T>
tuple<matrix<_T>,matrix<_T>> QR(const matrix<_T>& m){
	auto M=m;
	matrix<_T> v(0,0);
	int row(size(m).get(0)),col(size(m).get(1));
	auto H=matrix<_T>(eye(row));
	_T nrm;
	matrix<_T> submtr(0,0);
	for(int i=0;i<MIN(row,col);++i){
		v=M.get(linspace(i,row-1),linspace(i));
		v[0][0]+=SIGN(v[0][0])*vecnorm(v).get(0);
		nrm=pow(vecnorm(v).get(0),2);
		if(IS_NON_ZERO(nrm)){
			submtr=M.get(linspace(i,row-1),linspace(i,col-1));
			M.set(linspace(i,row-1),linspace(i,col-1),submtr-2/nrm*(v^(v.T()^submtr)));
			submtr=H.get(linspace(0,row-1),linspace(i,row-1));
			H.set(linspace(0,row-1),linspace(i,row-1),submtr-2/nrm*((submtr^v)^(v.T())));
		}
	}
	return make_tuple(H,M);
}

#endif



























