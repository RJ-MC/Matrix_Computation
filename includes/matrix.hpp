#ifndef matrix_cpp
#define matrix_cpp

#include "CONSTANTS.h"
#include "matrix.h"

#include <complex>
#include <cmath>
#include <string>
#include <initializer_list>
#include <type_traits>

using namespace std;

template <class _T>
matrix<_T>::matrix(int row, int col, initializer_list<_T>list){
    sze[0]=row;
    sze[1]=col;
	mtr = new _T*[row];
	for (int i = 0; i < row; ++i)
		mtr[i] = new _T[col]();
	int _row_ind = 0, _col_ind = -1;
	for (auto elem : list)
		if (++_col_ind < col)
			mtr[_row_ind][_col_ind] = elem;
   		else if (++_row_ind < row)
			mtr[_row_ind][_col_ind=0] = elem;
		else
			break;
}

template <class _T>
matrix<_T>::matrix(const matrix<int>& Size, initializer_list<_T>list){
    sze[0]=Size.get(0);
    sze[1]=Size.get(1);
	mtr = new _T*[sze[0]];
	for (int i = 0; i < sze[0]; ++i)
		mtr[i] = new _T[sze[1]]();
	int _row_ind = 0, _col_ind = -1;
	for (auto elem : list)
		if (++_col_ind < sze[1])
			mtr[_row_ind][_col_ind] = elem;
   		else if (++_row_ind < sze[0])
			mtr[_row_ind][_col_ind=0] = elem;
		else
			break;
}

template <class _T>
matrix<_T>::~matrix() {
	if (sze[0]) {
		for (int i = 0; i < sze[0]; ++i)
			delete[] mtr[i];
		delete[] mtr;
	}
}

template <class _T>
matrix<_T>::matrix(const matrix<_T> &M) {
	sze[0]=M.sze[0];
	sze[1]=M.sze[1];
	mtr = new _T*[sze[0]];
	for (int i = 0; i < sze[0]; ++i) {
		mtr[i] = new _T[sze[1]];
		for (int j = 0; j < sze[1]; ++j)
			mtr[i][j] = M.get(i, j);
	}
}
template <class _T>
template <class _U>
matrix<_T>::matrix(const matrix<_U> &M){
	sze[0]=size(M).get(0);
	sze[1]=size(M).get(1);
	mtr = new _T*[sze[0]];
	for (int i = 0; i < sze[0]; ++i) {
		mtr[i] = new _T[sze[1]];
		for (int j = 0; j < sze[1]; ++j)
			mtr[i][j] = M.get(i, j);
	}
}

template <class _T>
const _T matrix<_T>::get(int i, int j) const {
	return mtr[i][j];
}

template <class _T>
const _T matrix<_T>::get(int ind) const {
	return mtr[ind - ind / sze[0] * sze[0]][ind / sze[0]];
}

template <class _T>
ostream& operator << (ostream &os, const matrix<_T>& m) {
	os.setf(ios::fixed, ios::floatfield);
	os.precision(2);
	if (m.sze[0]*m.sze[1]) {
		for (int i = 0; i < m.sze[0]; ++i) {
			if(i==0)
				os << "©°";
			else if(i==m.sze[0]-1)
				os<< "©¸";
			else
				os<< "©¦";
			for (int j = 0; j < m.sze[1]; ++j)
				os << " " << setw(15) << right << m.get(i, j);
			if(i==0)
				os << "©´";
			else if(i==m.sze[0]-1)
				os<< "©¼";
			else
				os<< "©¦";
			os<<"\n";
		}
	}
	else
		os << "Empty-" << m.sze[0] << "x" << m.sze[1] << "-matrix\n";
	return os;
}

template <class _T>
matrix<_T> matrix<_T>::T(bool _conj) const {
	matrix<_T> m(sze[1], sze[0]);
	if(_conj)
		for (int i = 0; i < sze[1]; ++i)
			for (int j = 0; j < sze[0]; ++j)
				m[i][j] = conj(mtr[j][i]);
	else
		for (int i = 0; i < sze[1]; ++i)
			for (int j = 0; j < sze[0]; ++j)
				m[i][j] = mtr[j][i];
	return m;
}

template <class _T>
_T* matrix<_T>::operator [] (int _row) {
	return mtr[_row];
}

template <class _T>
matrix<int>& ind2sub(int ind, const matrix<_T>& size) {
	return matrix<int> (2, 1, {ind - ind / size.get(0,0) * size.get(0,0),ind / size.get(0,0)});
}

template <class _T>
int numel(const matrix<_T>& M) {
	return M.sze[0] * M.sze[1];
}

template <class _T>
matrix<int>& ind2sub(const matrix<int>& ind, const matrix<_T>& size) {
    int n(numel(ind));
    
	matrix<int>* res(new matrix<int> (2, n));
	for(int i=0;i<n;i++){
	    res->mtr[i][0]=ind / size.get(0,0);
	    res->mtr[i][1]=ind - ind / size.get(0,0) * size.get(0,0);
	}
	    
	return *res;
}

template <class _T>
matrix<int>& size(const matrix<_T>& M){
    return *(new matrix<int>(1,2,{M.sze[0],M.sze[1]}));
}

template <class _T>
matrix<_T>& matrix<_T>::operator[] (const matrix<int>& ind) const{
    matrix<_T>* M=new matrix<_T>(size(ind),{});
    for(int i=0;i<size(ind).get(0);++i)
        for(int j=0;j<size(ind).get(1);++j)
            (*M)[i][j]=this->get(ind.get(i,j));
	return *M;
}

template <class _T>
matrix<_T>& matrix<_T>::get (const matrix<int>& ind) const {
	return this->operator[](ind);
}

template <class _T>
matrix<_T>& matrix<_T>::get (const matrix<int>& row, const matrix<int>& col) const{
    matrix<_T>* M=new matrix<_T>(numel(row),numel(col));
    for(int i=0;i<numel(row);i++)
        for(int j=0;j<numel(col);j++)
            (*M)[i][j]=this->get(row.get(i),col.get(j));
	return *M;
}

template<class _U, class _V> bool operator == (const matrix<_U>& M1, const matrix<_V>& M2){
    if(M1.sze[0]!=M2.sze[0]||M1.sze[1]!=M2.sze[1])
        return 0;
    for(int i=0;i<M1.sze[0];i++)
        for(int j=0;j<M1.sze[1];j++)
            if(M1.mtr[i][j]!=M2.mtr[i][j])
                return 0;
    return 1;
}

template<class _U, class _V> bool operator != (const matrix<_U>& M1, const matrix<_V>& M2){
    return !(M1==M2);
}

template <class _T>
template <class _U>
matrix<_T>& matrix<_T>::set(const matrix<int>& ind, const matrix<_U>& val) {
	int _ind;
	if (numel(ind) == numel(val)) 
		for (int i = 0; i < numel(ind); ++i){
				_ind=ind.get(i);
				mtr[_ind-_ind/val.sze[0]*val.sze[0]][_ind/val.sze[0]] = _T(val.get(i));
		}
	else
		throw(string("[Error] set : Sizes do not match."));
	
	return *this;
}

template <class _T>
matrix<_T>& reshape(const matrix<_T>& M, const matrix<int>& size) {
    if(numel(size)<2)
        throw(string("[Error] reshape : Size matrix should have exactly 2 elements."));
	matrix<_T>*v=new matrix<_T>(size.get(0), size.get(1));

	int i(-1);

	while ((++i) < MIN(size.get(0)*size.get(1),numel(M)))
		(*v)[i-(i/size.get(0))*size.get(0)][i/size.get(0)] = M.get(i);

	return *v;
}

template <class _T>
matrix<_T>& matrix<_T>::operator=(const matrix<_T>& M){
    if(sze[0]!=M.sze[0]||sze[1]!=M.sze[1]){
        this->~matrix<_T>();
        new(this) matrix(M.sze[0],M.sze[1]);
    }
    for(int i=0;i<sze[0];i++)
        for(int j=0;j<sze[1];j++)
            mtr[i][j]=M.get(i,j);
    return *this;
}

template <class _T>
template <class _U>
matrix<_T>& matrix<_T>::operator=(const matrix<_U>& M){
    if(sze[0]!=size(M).get(0) || sze[1]!=size(M).get(1)){
        this->~matrix<_T>();
        new(this) matrix(size(M),{});
    }
    for(int i=0;i<sze[0];i++)
        for(int j=0;j<sze[1];j++)
            mtr[i][j]=M.get(i,j);
    return *this;
}

#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)

#define MATRIX_OPS(OP)\
template <class _T, class _U>\
auto operator OP(const matrix<_T>& m1, const matrix<_U>& m2)->matrix<typename remove_const<decltype(_T(1) OP _U(1))>::type>& {\
	if(size(m1)!=size(m2))\
		throw(string(STRINGIFY([Error] Operator OP : Size do not match.)));\
	auto *r= new\
		matrix<typename remove_const<decltype(_T(1) OP _U(1))>::type> (m1);\
	for (int i = 0; i < size(m1).get(0); ++i)\
		for (int j = 0; j < size(m1).get(1); ++j)\
			r->mtr[i][j] = r->mtr[i][j] OP m2.get(i,j);\
	return *r;\
}\
template <class _T, class _U>\
auto operator OP(const matrix<_T>& m1, const _U& m2)->matrix<typename remove_const<decltype(_T(1) OP _U(1))>::type>& {\
	auto  *r= new matrix<typename remove_const<decltype(_T(1) OP _U(1))>::type> (m1);\
	for (int i = 0; i < size(m1).get(0); ++i)\
		for (int j = 0; j < size(m1).get(1); ++j)\
			r->mtr[i][j] = r->mtr[i][j] OP m2;\
	return *r;\
}\
template <class _T, class _U>\
auto operator OP(const _T& m1, const matrix<_U>& m2)->matrix<typename remove_const<decltype(_T(1) OP _U(1))>::type>& {\
	auto *r= new matrix<typename remove_const<decltype(_T(1) OP _U(1))>::type> (m2);\
	for (int i = 0; i < size(m2).get(0); ++i)\
		for (int j = 0; j < size(m2).get(1); ++j)\
			r->mtr[i][j] = m1 OP r->mtr[i][j];\
	return *r;\
}
MATRIX_OPS(+)
MATRIX_OPS(-)
MATRIX_OPS(*)
MATRIX_OPS(/)
#undef MATRIX_OPS

template <class _T, class _U>
auto operator ^(const matrix<_T>& m1, const _U& m2)->matrix<typename remove_const<decltype(pow(_T(1),_U(1)))>::type>& {
	auto  *r= new matrix<typename remove_const<decltype(pow(_T(1),_U(1)))>::type>(m1);
	for (int i = 0; i < size(m1).get(0); ++i)
		for (int j = 0; j < size(m1).get(1); ++j)
			(*r)[i][j] = pow((*r)[i][j],m2);
	return *r;
}
template <class _T, class _U>
auto operator ^(const _T& m1, const matrix<_U>& m2)->matrix<typename remove_const<decltype(pow(_T(1),_U(1)))>::type>& {
	auto *r= new matrix<typename remove_const<decltype(pow(_T(1),_U(1)))>::type> (m2);
	for (int i = 0; i < size(m2).get(0); ++i)
		for (int j = 0; j < size(m2).get(1); ++j)
			(*r)[i][j] = pow(m1,(*r)[i][j]);
	return *r;
}

template <class _T, class _U>
auto operator ^(const matrix<_T>& m1, const matrix<_U>& m2)->matrix<typename remove_const<decltype(_T(1) * _U(1))>::type>& {
	if(m1.sze[1]!=m2.sze[0])
		throw(string("STRINGIFY([Error] Operator ^ : Size do not match."));
	auto* r= new 
		matrix<typename remove_const<decltype(_T(1) * _U(1))>::type>(size(m1).get(0),size(m2).get(1));
	_T s;
	for(int i=0;i<m1.sze[0];++i)
	    for(int k=0;k<m1.sze[1];++k){
	        s=m1.mtr[i][k];
	        for(int j=0;j<m2.sze[1];++j)
	            r->mtr[i][j]+=s*m2.mtr[k][j];
		}
	return *r;
}

#endif
