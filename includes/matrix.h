#ifndef matrix_h
#define matrix_h

#include <iostream>
#include <tuple>
#include <iomanip>
#include <initializer_list>
using namespace std;

template <class _T>
class matrix {
public:

	int size[2];
	_T** mtr;

	matrix(int row = 0, int col = 0, initializer_list<_T>list = initializer_list<_T>());
	matrix(const matrix<int>&, initializer_list<_T>list);
	matrix(const matrix<_T> &);
	~matrix();
	
	_T* operator[] (int _row);
	const _T get(int) const;
	const _T get(int, int) const;
	matrix<_T>& get(const matrix<int>&) const;
	matrix<_T>& get(const matrix<int>&,const matrix<int>&) const;
	matrix<_T>& operator[] (const matrix<int>&) const;
	
	template <class _U>
	matrix<_T>& set(const matrix<int>&, const matrix<_U>&);
	
	matrix T(bool _conj = true) const;
	
	template <class _U>
	operator matrix<_U>() const;
	
	matrix<_T>& operator=(const matrix<_T>& );
	template <class _U>
	matrix<_T>& operator=(const matrix<_U>& );
	
};

//What happens to me here?
#include "matrix.cpp"

#endif
