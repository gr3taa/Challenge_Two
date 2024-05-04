# Challenge_Two
Second challenge of the advanced programming for scientific computing. The task is to operate with large sparse matrices.

The sparse matrices can be stored in different ways. The techniques implemented in this project are:

+ **Compressed storage techniques:** 

this is the most efficient way to store the sparse matrix in terms of memory and of the computational efficiency of basic operations. t. But
they do not allow to change the pattern of sparsity.

+ **Uncompressed storage techniques:** 

it is possible to add a non zero elements easily.

Another distinction is given by the storage ordering, which can be either row-wise or column-wise

## Organization of the code

In the folder `src` there is:
+ [main.cpp](https://github.com/gr3taa/Challenge_Two/blob/main/src/main.cpp)

In the folder `include` there are:
+ [matrix_class.hpp](https://github.com/gr3taa/Challenge_Two/blob/main/include/matrix_class.hpp)
+ [function_implementation.hpp](https://github.com/gr3taa/Challenge_Two/blob/main/include/function_implementation.hpp)
+ [friend_functions.hpp](https://github.com/gr3taa/Challenge_Two/blob/main/include/friend_functions.hpp)
+ [Doxyfile](https://github.com/gr3taa/Challenge_Two/blob/main/include/Doxyfile)
+ [chrono.hpp](https://github.com/gr3taa/Challenge_Two/blob/main/include/chrono.hpp)

In the folder `file` there are:
+ [lnsp_131.mtx](https://github.com/gr3taa/Challenge_Two/blob/main/file/lnsp_131.mtx)
+ [input.txt](https://github.com/gr3taa/Challenge_Two/blob/main/file/input.txt)

In the main folder there are:
+ [Makefile](https://github.com/gr3taa/Challenge_Two/blob/main/Makefile)
+ [README.md](https://github.com/gr3taa/Challenge_Two/blob/main/README.md)


In `matrix.hpp` there is implemented the matrix class and the two enum classes used in the code. 

The matrix class is a template class, the template parameters are the type of the elements and an enumerator the indicates the storage ordering (row-wise or column-wise).

The enum classes are:
+ StorageOrder: it could be `Row_wise` or `Column_wise`.
+ Norm: it could be `One` or `Infinity` or `Frobenius`.


In `function_implementation.hpp` there are implemented the functions of matrix class.

The matrix_class functions are:
+ `resize()`
+ `compress()`
+ `uncompress()`
+ `is_compressed()`
+ `operator()`: there are two versions, const and non const version
+ `get_n_rows()`
+ `get_n_cols()`
+ `get_nnz()`
+ `set_nnz()`
+ `norm_one()`
+ `norm_infinity()`
+ `norm_frobenius()`
+ `norm()`
+ `print_matrix()`
+ `print_compressed_matrix()`
+  `read_matrix()`

In `friend_functions.hpp` there are implemented the friend functions of matrix class.
The friends functions are:
+ `operator<`: it have to be overloaded because if the storage order is `Column_wise` we need a column-major order fore the map.
+ `operator*`: one is for the matrix-vector product and the other for matrix-matrix product.

The description of all these functions can be found in the doxygen documentation running `make doc`.

`chrono.hpp` allows to time the operations.

`lnsp_131.mtx` is the file that contains the sparse matrix.

`input.txt` in this file there are written the storage order and the type of the norm computed.

## Running locally

Run `make` to compile the program and then `./program`. 

Run `make doc` to generate the documentation of the project.

The user can choose the storage order and the norm to be calculated modifying the `input.txt`.



