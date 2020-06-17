#ifndef matrix_h
#define matrix_h

#include <iostream>
#include <cmath>
#include <iomanip>
#include <initializer_list>
using namespace std;


template <class _T>
class matrix {

	int sze[2];
	_T** mtr;
	
public:

	matrix(int row, int col, initializer_list<_T>list = initializer_list<_T>());
	matrix(const matrix<int>&, initializer_list<_T>list);
	matrix(const matrix<_T> &);
	template <class _U>
	matrix(const matrix<_U> &);
	~matrix();
	
	_T* operator[] (int _row);
	const _T get(int) const;
	const _T get(int, int) const;
	matrix<_T>& get(const matrix<int>&) const;
	matrix<_T>& get(const matrix<int>&,const matrix<int>&) const;
	matrix<_T>& operator[] (const matrix<int>&) const;
	
	template <class _U>
	matrix<_T>& set(const matrix<int>&, const matrix<_U>&);
	template <class _U>
	matrix<_T>& set(const matrix<int>&, const matrix<int>&, const matrix<_U>&);
	
	matrix T(bool _conj = true) const;
	
	matrix<_T>& operator=(const matrix<_T>& );
	template <class _U>
	matrix<_T>& operator=(const matrix<_U>& );
	
	template<class _U>
	friend ostream& operator << (ostream &, const matrix<_U>&);
	
	template<class _U>
	friend istream& operator >> (istream& input, matrix<_U>&);
	
	template<class _U, class _V>
	friend bool operator == (const matrix<_U>&, const matrix<_V>&);
	
	template <class _U>
	friend int numel(const matrix<_U>&);
	
	template <class _U>
	friend matrix<int>& size(const matrix<_U>& );
	
	#define MATRIX_OPS(OP)\
	template <class _V, class _U>\
	friend auto operator OP(const matrix<_V>& m1, const matrix<_U>& m2)->matrix<typename remove_const<decltype(_V(1) OP _U(1))>::type>&;\
	template <class _V, class _U>\
	friend auto operator OP(const matrix<_V>& m1, const _U& m2)->matrix<typename remove_const<decltype(_V(1) OP _U(1))>::type>&;\
	template <class _V, class _U>\
	friend auto operator OP(const _V& m1, const matrix<_U>& m2)->matrix<typename remove_const<decltype(_V(1) OP _U(1))>::type>&;
	MATRIX_OPS(+)
	MATRIX_OPS(-)
	MATRIX_OPS(*)
	MATRIX_OPS(/)
	#undef MATRIX_OPS
	
	template <class _V, class _U>
	friend auto operator ^(const matrix<_V>& m1,const  matrix<_U>& m2)->matrix<typename remove_const<decltype(_V(1) * _U(1))>::type>&;
	template <class _V, class _U>
	friend auto operator ^(const matrix<_V>& m1, const _U& m2)->matrix<typename remove_const<decltype(pow(_V(1),_U(1)))>::type>&;
	template <class _V, class _U>
	friend auto operator ^(const _V& m1, const matrix<_U>& m2)->matrix<typename remove_const<decltype(pow(_V(1),_U(1)))>::type>&;
	
	#define MATRIX_REL_OPS(OP)\
	template <class _V, class _U>\
	friend matrix<bool>& operator OP(const matrix<_V>& m1, const _U& m2);\
	template <class _V, class _U>\
	friend matrix<bool>& operator OP(const _V& m1, const matrix<_U>& m2);
	MATRIX_REL_OPS(<)
	MATRIX_REL_OPS(<=)
	MATRIX_REL_OPS(>)
	MATRIX_REL_OPS(>=)
	MATRIX_REL_OPS(==)
	#undef MATRIX_REL_OPS
	
	friend matrix<int>& eye(int);
	friend matrix<int>& ones(const matrix<int> &);
	
};

//What happens to me here?
#include "matrix.hpp"
#include "spc_mtr.hpp"
#include "mtr_fcn.hpp"

#endif
