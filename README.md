# Challenge_Two
Second challenge of the advanced programming for scientific computing. The task is to operate with large matrices.

## Organization of the code

In the folder `src` there are:
+ [main.cpp](https://github.com/gr3taa/Challenge_Two/blob/main/src/main.cpp)
+ [matrix.hpp](https://github.com/gr3taa/Challenge_Two/blob/main/src/matrix.hpp)
+ [function_implementation.hpp](https://github.com/gr3taa/Challenge_Two/blob/main/src/function_implementation.hpp)
+ [chrono.hpp](https://github.com/gr3taa/Challenge_Two/blob/main/src/chrono.hpp)
+ [Insp_131.mtx](https://github.com/gr3taa/Challenge_Two/blob/main/lnsp_131.mtx)
+ [Doxyfile](https://github.com/gr3taa/Challenge_Two/blob/main/Doxyfile)

In `matrix.hpp` there is implemented the matrix class. It is a template class, the template parameters are the type of the elements and an enumerator the indicates the storage ordering (row-wise or column-wise).
The sparse matrix can be stored in two different storage techniques:
+ **Compressed storage techniques:** 
this is the most efficient way to store the sparse matrix in terms of memory and of the computational efficiency of basic operations. The `compress()` function converts the internal storage from the map to a compresses sparse row (or columns, depending on the matrix ordering). You can check if a matrix is compressed through `is_compressed()`.
+ **Uncompressed storage techniques:** 
it is possible to add a non zero elements through the call operator.

In `function_implementation.hpp` there are implemented the functions of matrix class and its friend function:
+ product function between a matrix and a vector.
+ comparison operator.


`chrono.hpp` allows to time the operations.

`lnsp_131.mtx` is the file that contains the sparse matrix.

## Running locally

Run `make` to compile the program and then `./program`. 
Run `make doc` to generate the documentation of the project.



