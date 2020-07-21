# Matrix_Computations
Implementation of *Matrix Computations (4th)* via `C++`.

## Class `matrix`

<details>
<summary>Usage</summary>

```cpp
#include "includes/matrix.h"
```
</details>
<details>
<summary>Basic Opreations</summary>

-   <details>
    <summary>Create matrix</summary>

    ```cpp
    matrix<_T> M(n,m);          // creates M as a n-by-m matrix of type _T with 0 ...
    matrix<_T> M(n,m,{...});    // or elements listed in {...} with column-first order
    matrix<_T> M(Size,{...});   // Size is a matrix of at least 2 element.
    ```
    - Examples
        ```cpp
        matrix<int> M(2,2);
        matrix<int> M(2,2,{1,2,3});
        matrix<int> M(2,2,{1,2,3,4,5});
        matrix<int> M(matrix<int>(2,1,{2,2}),{});
        ```
        define
        <p align="center"><img src="https://latex.codecogs.com/png.latex?\begin{bmatrix}0&0\\0&0\end{bmatrix},\quad\begin{bmatrix}1&2\\3&0\end{bmatrix},\quad\begin{bmatrix}1&2\\3&4\end{bmatrix},\quad\begin{bmatrix}0&0\\0&0\end{bmatrix}."></p>
    </details>



-   <details>
    <summary>Get/set single element</summary>

    There are two ways to get element of a matrix, via **index** or **subscripts**. Index starts from 0 at upper-left corner, and increases 1 by moving down or moving to the top of the next column. Subscripts starts from (0,0) at upper-left corner.
    ```cpp
    M.get(i);        // get via index i.
    M.get(i,j);      // get via subscripts (i,j).
    M[i][j]=...;     // set via subscrpits (i,j).
    ```
    </details>

-   <details>
    <summary>Get submatrix</summary>

    ```cpp
    M[I];           // I is an index matrix.
    M.get(I);
    M.get(R,L);     // R is index matrix of row and L is of column.
    ```
    - Example
        ```cpp
        matrix<int>(2,2,{1,2,3,4}).submtr(matrix<int>(2,1,{0,0}),matrix<int>(1,2,{1,0}));   // [2,1;2,1]
        ```
    </details>

-   <details>
    <summary>Set submatrix</summary>

    ```cpp
    M.set(I,V);         // I is index matrix, V is value matrix with same element number of I.
    ```
    </details>
</details>

<details>
<summary>Advanced Opreations</summary>

```cpp
M.T();              // return conjugate of M.
M.T(0);             // return transpose of M.
M*N;                // point-wise multiplication of two same-size matricies.
M/N;                // point-wise division of two same-size matricies.
M^N;                // matrix multiplication.
```
</details>

<details>
<summary>Special matricies</summary>

```cpp
eye(n);             // create an identity matrix of size n.
ones(Size);         // create an all-1 matrix of size matrix Size.
zeros(Size);        // create an all-0 matrix of size matrix Size.
linspace(n,m,d);    // create a row-1 matrix, starting from n,
                    // increasing by d, and ending with the number
                    // whose next step will be greater than m.
```
</details>

<details>
<summary>Matrix functions</summary>

```cpp
min(M);             // return a row-1 matrix consisting of minimum value of each column.
max(M);             // return a row-1 matrix consisting of maximum value of each column.
chs(M,f);           // return a row-1 matrix consisting of maximum value of each column
                    // by the comparing function f.
                    // Example: For int matrix M,
                    //      chs(M,(bool(*)(const int&,const int&))
                    //      ([](const int& a,const int& b)->bool{return a>=b;}));)=max(M).
size(M);            // return a 2-by-1 matrix indicating row and column numbers.
numel(M);           // return the number of elements.
reshape(M,Size);    // Size is an at-least-2-element matrix.
vecnorm(M,p=2);     // return a row-1 matrix consisting of lp norm of each column.
diag(M,n=0);        // return the n-offset diagonal of M as a row-1 matrix.
utri(M,n=0);
ltri(M,n=0);        // return the upper/lower part of M with offset n;
apply(M,f);         // apply f to each element of M.
find(M);            // for a bool M, return a row-1 matrix consisting of the location of all trues.
find(M,false);      // for a bool M, return a row-1 matrix consisting of the location of the first true.
diff(M);            // return difference of each column: M(2:end,:)-M(1:end-1,:).
det(M);             // return the determinant.
abs(M);             // return abs of M.
cmb(M1,M2,_horizontal=true)
                    // combine M1 and M2 in the indicated order.
```
</details>

<details>
<summary>Miscellanies</summary>

```cpp
simplex(M);         // solve linear programming problem via two-phase method.
```
</details>

## `Chapter 1` &emsp; Matrix Multiplication

<details>
<summary>Matrix functions</summary>

```cpp
fft2(M);        // conduct FFT to each column of a power-of-2-row M.
fht(M);         // conduct FHT to each column of a power-of-2-row M.
```
</details>

## `Chapter 5` &emsp; Orthogonalization and Least Squares

<details>
<summary>Matrix functions</summary>

```cpp
QR(M);          // Solve QR decomposition of M.
                // Q = get<0>(QR(M)); R = get<1>(QR(M));
```
</details>

## `Chapter 7` &emsp; Unsymmetric Eigenvalue Problems

<details>
<summary>Matrix functions</summary>

```cpp
eigval(M);      // return all eigenvalues of M.
roots(M);       // return zeros of polynomial of coefficients M.
```
</details>