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
	matrix(const matrix<_T> &);
	~matrix();
	
	const _T get(int) const;
	const _T get(int, int) const;
	
	matrix T(bool _conj = true) const;
	_T* operator[] (int _row);
	
	matrix<_T>& operator[] (const matrix<int>&);
	
	int numel() const;
	
	/*;
	matrix utri(int = 0);
	matrix ltri(int = 0);
	*/
	
	matrix<_T>& reshape(const matrix<int>&) const;
	
	template <class _U>
	operator matrix<_U>() const;
	
	matrix<_T>& operator=(const matrix<_T>& );
	
	template <class _U>
	matrix<_T>& operator=(const matrix<_U>& );
	
	template <class _U>
	void set(const matrix<int>&, const matrix<_U>&);
	
	template<class _U, class _V>
	friend bool operator == (const matrix<_U>&, const matrix<_V>&);
	template<class _U, class _V>
	friend bool operator != (const matrix<_U>&, const matrix<_V>&);
	
	template <class _U>
	friend ostream& operator << (ostream &, const matrix<_U>&);

	matrix<_T>& submtr(const matrix<int>&, const matrix<int>&);

};

//What happens to me here?
#include "matrix.cpp"

#endif

