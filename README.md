# Matrix_Computation
Implementation of *Matrix Computations (4th)* by `C++`.

## Class `matrix`

### Usage
```cpp
#include "includes/matrix.h"
```

### Basic Operations
#### Create matrix

- Syntax
	```cpp
	matrix<_T> M(n,m);
	matrix<_T> M(n,m,{...});
	```
	creates `M` as a `n`-by-`m` matrix of type `_T` with `0`, or elements listed in `{...}` with column-first order.
	
- Examples
	```cpp
	matrix<int> M(2,2);
	```
	defines matrix
	<p align="center"><img src="https://latex.codecogs.com/png.latex?\begin{bmatrix}0&0\\0&0\end{bmatrix},"></p>
	and 
	
	```cpp
	matrix<int> M(2,2,{1,2,3});
	matrix<int> M(2,2),{1,2,3,4,5});
	```
	define
	<p align="center"><img src="https://latex.codecogs.com/png.latex?\begin{bmatrix}1&2\\3&0\end{bmatrix},\quad\begin{bmatrix}1&2\\3&4\end{bmatrix}."></p>

#### Get/set single element

There are two ways to get element of a matrix, via **index** or **subscripts**. Index starts from 0 at upper-left corner, and increases 1 by moving down or moving to the top of the next column. Subscripts starts from (0,0) at upper-left corner.
	
- Syntax
	```cpp
	M.get(i);        // get via index i.
	M.get(i,j);      // get via subscripts (i,j).
	M[i][j]=...;     // set via subscrpits (i,j).
	```
#### Get submatrix
- Syntax
	```cpp
	M[I];           // I is an index matrix.
	M.get(I);
	M.get(R,L);     // R is index matrix of row and L is of column.
	                // Numbers of elements of R and L should coincide.
	```
- Example
	```cpp
	matrix<int>(2,2,{1,2,3,4}).submtr(matrix<int>(2,1,{0,0}),matrix<int>(1,2,{1,0}));	// [2,1;2,1]
	```
#### Set submatrix
- Syntax
	```cpp
	M.set(I,V);	// I is index matrix, V is value matrix with same element number of I.
	```
#### Get size
- Syntax
	```cpp
	size(M);  	// return a 2-by-1 matrix indicating row and column numbers.
	```
#### Get number of elements
- Syntax
	```cpp
	numel(M);	// return the number of elements.
	```
### Advanced Operations

#### Transpose/Conjugate
- Syntax
	```cpp
	M.T();   	// return conjugate of M.
	M.T(0);   	// return transpose of M.
	```
#### Reshape
- Syntax
	```cpp
	reshape(M,Size);	// Size is a 2-element matrix.
	```
