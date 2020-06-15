#ifndef matrix_cpp
#define matrix_cpp

#include "CONSTANTS.h"
#include "matrix.h"

#include <cmath>
#include <string>
#include <initializer_list>
#include <type_traits>

#include <complex>
#include <algorithm>
using namespace std;

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
matrix<_T>::matrix(int row, int col, initializer_list<_T>list){
    size[0]=row;
    size[1]=col;
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
    size[0]=Size.get(0);
    size[1]=Size.get(1);
	mtr = new _T*[size[0]];
	for (int i = 0; i < size[0]; ++i)
		mtr[i] = new _T[size[1]]();
	int _row_ind = 0, _col_ind = -1;
	for (auto elem : list)
		if (++_col_ind < size[1])
			mtr[_row_ind][_col_ind] = elem;
   		else if (++_row_ind < size[0])
			mtr[_row_ind][_col_ind=0] = elem;
		else
			break;
}

template <class _T>
matrix<_T>::~matrix() {
	if (size[0]) {
		for (int i = 0; i < size[0]; ++i)
			delete[] mtr[i];
		delete[] mtr;
	}
}

template <class _T>
matrix<_T>::matrix(const matrix<_T> &M) {
	size[0]=M.size[0];
	size[1]=M.size[1];
	mtr = new _T*[size[0]];
	for (int i = 0; i < size[0]; ++i) {
		mtr[i] = new _T[size[1]];
		for (int j = 0; j < size[1]; ++j)
			mtr[i][j] = M.get(i, j);
	}
}

template <class _T>
template <class _U>
matrix<_T>::matrix(const matrix<_U> &M) {
	size[0]=M.size[0];
	size[1]=M.size[1];
	mtr = new _T*[size[0]];
	for (int i = 0; i < size[0]; ++i) {
		mtr[i] = new _T[size[1]];
		for (int j = 0; j < size[1]; ++j)
			mtr[i][j] = M.get(i, j);
	}
}

template <class _T>
const _T matrix<_T>::get(int i, int j) const {
	return mtr[i][j];
}

template <class _T>
const _T matrix<_T>::get(int ind) const {
	return mtr[ind - ind / size[0] * size[0]][ind / size[0]];
}

template <class _T>
ostream& operator << (ostream &os, const matrix<_T>& m) {
	os.setf(ios::fixed, ios::floatfield);
	os.precision(2);
	if (m.size[0]*m.size[1]) {
		for (int i = 0; i < m.size[0]; ++i) {
			if(i==0)
				os << "©°";
			else if(i==m.size[0]-1)
				os<< "©¸";
			else
				os<< "©¦";
			for (int j = 0; j < m.size[1]; ++j)
				os << " " << setw(15) << right << m.get(i, j);
			if(i==0)
				os << "©´";
			else if(i==m.size[0]-1)
				os<< "©¼";
			else
				os<< "©¦";
			os<<"\n";
		}
	}
	else
		os << "Empty-" << m.size[0] << "x" << m.size[1] << "-matrix\n";
	return os;
}

template <class _T>
matrix<_T> matrix<_T>::T(bool _conj) const {
	matrix<_T> m(size[1], size[0]);
	if(_conj)
		for (int i = 0; i < size[1]; ++i)
			for (int j = 0; j < size[0]; ++j)
				m[i][j] = conj(mtr[j][i]);
	else
		for (int i = 0; i < size[1]; ++i)
			for (int j = 0; j < size[0]; ++j)
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
	return M.size[0] * M.size[1];
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
    return *(new matrix<int>(1,2,{M.size[0],M.size[1]}));
}

template <class _T>
matrix<_T>& matrix<_T>::operator[] (const matrix<int>& ind) const{
    matrix<_T>* M=new matrix<_T>(ind.size[0],ind.size[1]);
    for(int i=0;i<ind.size[0];i++)
        for(int j=0;j<ind.size[1];j++)
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
    if(M1.size[0]!=M2.size[0]||M1.size[1]!=M2.size[1])
        return 0;
    for(int i=0;i<M1.size[0];i++)
        for(int j=0;j<M1.size[1];j++)
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
	if (numel(ind) && numel(val)) 
		for (int i = 0; i < numel(ind); ++i)
			for (int j = 0; j < val.numel(); ++j)
				mtr[i-i/val.size[0]*val.size[0]][i/val.size[0]] = _T(val.get(j));
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
    if(size[0]!=M.size[0]||size[1]!=M.size[1]){
        this->~matrix<_T>();
        new(this) matrix(M.size[0],M.size[1]);
    }
    for(int i=0;i<size[0];i++)
        for(int j=0;j<size[1];j++)
            mtr[i][j]=M.get(i,j);
    return *this;
}

template <class _T>
template <class _U>
matrix<_T>& matrix<_T>::operator=(const matrix<_U>& M){
    if(size[0]!=M.size[0]||size[1]!=M.size[1]){
        this->~matrix<_T>();
        new(this) matrix(M.size[0],M.size[1]);
    }
    for(int i=0;i<size[0];i++)
        for(int j=0;j<size[1];j++)
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
auto operator ^(matrix<_T>& m1, matrix<_U>& m2)->matrix<typename remove_const<decltype(m1.type * m2.type)>::type>& {
	if(m1.size[1]!=m2.size[0])
		throw(string("STRINGIFY([Error] Operator ^ : Size do not match."));
	matrix<typename remove_const<decltype(m1.type * m2.type)>::type>* r= new 
		matrix<typename remove_const<decltype(m1.type * m2.type)>::type>(m1.size[0],m2.size[1]);
	_T s;
	for(int i=0;i<m1.size[0];++i)
	    for(int k=0;k<m2.size[0];++k){
	        s=m1[i][k];
	        for(int j=0;j<m2.size[1];++j)
	            r->mtr[i][j]+=s*m2[k][j];
		}
	return *r;
}

/*
template <class _T>
tuple<matrix<_T>, matrix<_T>> qr(const matrix<_T>& m, const int& _start_row = 0) {
	matrix<_T> q(Id(m.row, _T(1))), r(m), H(q), x;
	matrix<int> ind_row, ind_q_col(seq(0, m.row - 1)), ind_r_col(seq(0, m.col - 1));
	for (int i = 0; i < (m.row > m.col ? m.col : m.row - 1) - _start_row; ++i) {
		ind_row = seq(i + _start_row, m.row - 1);
		x = r.submtr(ind_row, seq(i, i));
		if (norm(x) > EPS) {
			x[0][0] += norm(x)*SIGN(x[0][0]);
			H = Id(x.row) - 2 * x*x.T() / (x.T()*x);
			q.set(ind_row, ind_q_col, H*q.get(ind_row, ind_q_col));
			r.set(ind_row, ind_r_col, H*r.get(ind_row, ind_r_col));
		}
	}
	return make_tuple(q.T(), r.utri(-_start_row));
}

template <class _T, class _U>
auto operator / (const matrix<_T>& m1, const matrix<_U>& m2)->matrix<decltype(_U(0) / _T(1))> {
	matrix<decltype(_T(0)*_U(0))> r;
	if (m1.col == m2.col) {
		r = m1;
		auto qr_m2 = qr(m2);
		auto Q(get<0>(qr_m2)), R(get<1>(qr_m2));
		return r;
	}
	else if (m2.row == 1 && m2.col == 1) {
		_U den(m2.get(0, 0));
		r = 1 / den * m1;
	}
	else
		throw("[Error] Operator / : Sizes do not match.");
	return r;
}

template <class _T>
matrix<_T> rep(const matrix<_T>& _matrix, int _row, int _col) {
	matrix<_T> _rep(_matrix.row*_row, _matrix.col*_col);
	matrix<int> _row_seq(seq(0, _matrix.row - 1)), _col_seq;
	for (int i = 0; i < _row; ++i) {
		_col_seq = seq(0, _matrix.col - 1);
		for (int j = 0; j < _col; ++j) {
			_rep.set(_row_seq, _col_seq, _matrix);
			_col_seq = _col_seq + _matrix.col;
		}
		_row_seq = _row_seq + _matrix.row;
	}
	return _rep;
}

template <class _T>
matrix<int> sort(const matrix<_T>& m, bool increasing = 1) {
	matrix<int> ord(rep(seq(0, m.col - 1), m.row, 1));
	for (int i = 0; i < m.row; ++i)
		sort(&(ord[i][0]), &(ord[i][m.col]), [&ord, &m, i, increasing](int i1, int i2) {return increasing ^ (m.get(i, i1) > m.get(i, i2)); });
	return ord;
}

template <class _T>
matrix<int> max(const matrix<_T>& m) {
	matrix<int> ord(seq(0, m.col - 1)), ord_temp(seq(0, m.row - 1));
	for (int i = 0; i < m.col; ++i)
		ord[0][i] = *max_element(&(ord_temp[0][0]), &(ord_temp[0][m.col]), [&m, i](int i1, int i2) {return m.get(i1, i) < m.get(i2, i); });
	return ord;
}

template <class _T>
tuple<matrix<_T>, matrix<_T>> Hessen(const matrix<_T> &m) {
	// Return [Q,R] s.t. m=Q*R*Q'
	matrix<_T> q(Id(m.row, _T(1))), r(m), H(q), x(0, 0);
	if (m.row != m.col) {
		throw("[Error] Hessen : Not a square matrix.");
		return make_tuple(q, r);
	}
	matrix<int> ind_row(0, 0), ind_q_col(seq(0, m.col - 1)), ind_r_col(seq(0, m.col - 1));
	for (int i = 0; i < (m.row > m.col ? m.col : m.row - 1) - 1; ++i) {
		ind_row = seq(i + 1, m.row - 1);
		x = r.submtr(ind_row, seq(i, i));
		if (norm(x) > EPS) {
			x[0][0] += norm(x)*SIGN(x[0][0]);
			H = Id(x.row) - 2 / pow(norm(x), 2)*x*x.T();
			q.set(ind_row, ind_q_col, H*q.get(ind_row, ind_q_col));
			r.set(ind_row, ind_r_col, H*r.get(ind_row, ind_r_col));
			r.set(ind_r_col, ind_row, r.get(ind_r_col, ind_row)*H);
		}
	}
	return make_tuple(q.T(), r.utri(-1));
}

template <class _T>
matrix<complex<double>> fft(const matrix<_T> &y) {
	matrix<complex<double>> ffty(y.row,y.col);
	if (y.row)
		if((y.row & (y.row-1))==0) 
			if(y.row==1)
				return y.get(seq(0,0),seq(0,y.col-1));
			else{
				matrix<complex<double>> yt(fft(y.get(seq(0,y.row-1,2),seq(0,y.col-1)))),
										yb(fft(y.get(seq(1,y.row-1,2),seq(0,y.col-1))));
				matrix<complex<double>> z(dot_mult(rep(dot_pow(rep(
											matrix<complex<double>>(1,1,{complex<double>(cos(-2*PI/double(y.row)),sin(-2*PI/double(y.row)))})
											,1,y.row>>1),seq(0,(y.row>>1)-1)).T(false),1,y.col),yb));
				ffty.set(seq(0,(y.row>>1)-1),seq(0,y.col-1),yt+z);
				ffty.set(seq(y.row>>1,y.row-1),seq(0,y.col-1),yt-z);
			}
		else{
			int next_pow_2(1);
			while(next_pow_2<y.row)
				next_pow_2<<=1;
			matrix<complex<double>> pow_2_y(next_pow_2,y.col);
			pow_2_y.set(seq(0,y.row-1),seq(0,y.col-1),y);
			ffty=fft(pow_2_y).get(seq(0,y.row-1),seq(0,y.col-1));
		}
	return ffty;
}

template <class _T>
matrix<_T> fht(const matrix<_T> &y) {
	matrix<_T> fhty(y.row,y.col);
	if (y.row)
		if((y.row & (y.row-1))==0)
			if(y.row==1)
				return y.get(seq(0,0),seq(0,y.col-1));
			else{
				matrix<_T> z(fht(y.get(seq(0,(y.row>>1)-1),seq(0,y.col-1))));
				fhty.set(seq(0,y.row-1,2),seq(0,y.col-1),z+y.get(seq(y.row>>1,y.row-1),seq(0,y.col-1)));
				fhty.set(seq(1,y.row-1,2),seq(0,y.col-1),z-y.get(seq(y.row>>1,y.row-1),seq(0,y.col-1)));
			}
		else
			throw("[Error] fht : Only defined for matrix has 2-powered rows.");
	return fhty;
}*/

#endif
