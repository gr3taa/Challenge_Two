# Challenge_Two
Second challenge of the advanced programming for scientific computing. The task is to operate with large matrices.

## Organization of the code

In the folder `src` there are:
+ [main.cpp](https://github.com/gr3taa/Challenge_Two/blob/main/src/main.cpp)
+ [helper.hpp](https://github.com/gr3taa/Challenge_Two/blob/main/src/helper.hpp)
+ [matrix.hpp](https://github.com/gr3taa/Challenge_Two/blob/main/src/matrix.hpp)
+ [matrix.cpp](https://github.com/gr3taa/Challenge_Two/blob/main/src/matrix.cpp)
+ [chrono.hpp](https://github.com/gr3taa/Challenge_Two/blob/main/src/chrono.hpp)
+ [Insp_131.mtx](https://github.com/gr3taa/Challenge_Two/blob/main/src/Insp_131.mtx)

In `helper.hpp` there are implemented the friend function of matrix class and the norm method:
+ product function between a matrix and a vector
+ read_matrix that read the matrix written in matrix market format.
+ norm function that compute the one, infinity and Frobenius norm through the others function `norm_one`, `norm_infinity`,`norm_frobenius`.

In `matrix.hpp` and `matrix.cpp` there is implemented the matrix class. It is a template class, the template parameters are the type of the elements and an enumerator the indicates the storage ordering (row-wise or column-wise).
The sparse matrix can be stored in two different storage techniques:
+ **Compressed storage techniques:** 
this is the most efficient way to store the sparse matrix in terms of memory and of the computational efficiency of basic operations. The `compress()` function converts the internal storage from the map to a compresses sparse row (or columns, depending on the matrix ordering). You can check if a matrix is compressed through `is_compressed()`.
+ **Uncompressed storage techniques:** 
it is possible to add a non zero elements through the call operator.

`chrono.hpp` allows to time the matrix-vector product.

`Insp_131.mtx` is the file that contains the sparse matrix.


