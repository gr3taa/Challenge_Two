#ifndef matrix_class_hpp
#define matrix_class_hpp


#include<iostream>
#include<string>
#include<vector>
#include<array>
#include<map>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include<numeric>
#include<concepts>
#include <complex>
#include<type_traits>

using namespace std;

namespace algebra{
    
    /**
     * @brief Enum defining the storage order of the matrix.
     */
    enum class StorageOrder{
        Row_wise,
        Column_wise
    };

    /**
     * @brief Enum defining the type of norm.
     */
    enum class Norm{
        One,
        Infinity,
        Frobenius
    };

    using Index = array<size_t,2>;

    /**
     * @brief It is a dynamic matrix class template.
     * @tparam T is the type of elements.
     * @tparam ord is the storage order of the matrix, it can be Row_wise or Column_wise.
     */

    template<typename T, StorageOrder ord>
    class Matrix{
        private:
            // data is a map of indexes and elements of the matrix in case the matrix is not compressed
            map<Index, T> data;
            // I is the inner vector of indexes, II is the outer one in case the matrix is compressed
            vector<size_t> I, II; 
            // V is the vector of the values of nnz elements in case the matrix is compressed
            vector<T> V; 
            // n_rows is the number of rows of the matrix, n_cols is the number of columns, nnz_elem is the number of non zero elements
            size_t n_rows, n_cols, nnz_elem=0;
            // compressed is a boolean that is true if the matrix is compressed, false otherwise
            bool compressed; 
            /**
            * @brief Check if i and j are whitin the range of the matrix.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @param i is the index of the row.
            * @param j is the index of the column.
            * @return true if i and j are whitin the range of the matrix, false otherwise.
            */
            bool is_inbound(size_t i, size_t j) const;
            /**
            * @brief This function calculates the number of non zero elements of the matrix.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @return The number of non zero elements of the matrix.
            */
            size_t calculate_nnz() const;

        public:
            /**
             * @brief Constructor of the matrix class.
             * @param r is the number of rows of the matrix.
             * @param c is the number of columns of the matrix.
             */
            Matrix(size_t r, size_t c): n_rows(r), n_cols(c), compressed(false){
                nnz_elem = calculate_nnz();
                resize(r,c);

            };
            /**
            * @brief resizing the matrix.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @param i is the new size of rows.
            * @param j is the new size of colomns.
            */
            void resize(size_t i, size_t j);
            /**
            * @brief Converting the internal storage from a map to a compresses sparse row on comlumns depending on order.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            */
            void compress();
            /**
            * @brief Bringing back the matrix to the uncompressed state.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            */
            void uncompress(); 
            /**
            * @brief This function checks if the matrix is compressed or not.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @return true if the matrix is compressed.
            */
            bool is_compressed() const;
            /**
            * @brief Overloading call operator
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @param i is the index of the row.
            * @param j is the index of the column.
            * @return 0 if the element is within the range of the matrix but is not present.
            */
            T operator()(size_t i, size_t j) const; 
            /**
            * @brief Overloading call operator.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @param i is the index of the row.
            * @param j is the index of the column.
            * @return 
            */
            T& operator()(size_t i, size_t j); 
            /**
            * @brief Get the number of rows of the matrix.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @return the number of rows of the matrix.
            */  
            size_t get_n_rows() const;
            /**
            * @brief Get the number of columns of the matrix.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @return the number of columns of the matrix.
            */ 
            size_t get_n_cols() const;
            /**
            * @brief Get the number of non zero elements of the matrix.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @return the number of non zero elements of the matrix.
            */  
            size_t get_nnz() const;
            /**
            * @brief Set the number of non zero elements of the matrix.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @param nnz is the number of non zero elements.
            */ 
            void set_nnz(size_t nnz);
            /**
            * @brief this function computes the one norm.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @param mat the matrix of which we want to calculate the norm.
            * @return the one norm.
            */
            T norm_one() const;
            /**
            * @brief this function computes the infinity norm.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @param mat the matrix of which we want to calculate the norm.
            * @return the infinity norm.
            */
            T norm_infinity() const;
            /**
            * @brief this function computes the frobenius norm.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @param mat the matrix of which we want to calculate the norm.
            * @return the frobenius norm.
            */
            T norm_frobenius() const;
            
            /**
            * @brief this function calculates the norm of the matrix.
            * @tparam U is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @tparam n is the type of norm.
            * @param mat is the matrix that we want to calculate the norm.
            * @return the norm of the matrix.
            */
            template<Norm n>
            T norm() const;
            /**
             * @brief this function prints the matrix.
             * @tparam T is the type of elements.
             * @tparam order is the storage order of the matrix.
             */
            void print_matrix() const;
            /**
             * @brief this function prints the compressed matrix.
             * @tparam T is the type of elements.
             * @tparam order is the storage order of the matrix.
            */
            void print_compressed_matrix() const;

            /**
            * @brief this function reads the matrix written in matrix market format.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @param filename is the name of the file containing the matrix.
            * @return the matrix.
            */
            void read_matrix(const string & filename);

            /**
            * @brief Comparison operator according to the storage order.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @param lhs is the first index.
            * @param rhs is the second index.
            * @return true if lhs is less than rhs, false otherwise.
            */
            friend bool operator<(const Index & lhs, const Index & rhs);       

            /**
            * @brief this function does the matrix-vector product.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @param mat is the matrix that we want to multiply.
            * @param vec is the vector that we want to multiply.
            * @return the result of the matrix-vector product.
            */
            template<typename U, StorageOrder order>
            friend auto operator*(const Matrix<U,order> & mat, const vector<U> & vec);     

            /**
            * @brief this function does the matrix-matrix  product.
            * @tparam T is the type of elements.
            * @tparam order is the storage order of the matrix.
            * @param mat is the matrix that we want to multiply.
            * @param mat1 is the matrix that we want to multiply.
            * @return the result of the matrix-matrix product.
            */
            template<typename U, StorageOrder order>
            friend Matrix<U,order> operator*(const Matrix<U,order> & matl, const Matrix<U,order> & matr);     

    };
}

#include "function_implementation.hpp"
#include "friend_functions.hpp"
#endif 