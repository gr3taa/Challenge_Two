#ifndef matrix_hpp
#define matrix_hpp


#include<iostream>
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
     * @brief Matrix a dynamic matrix class template.
     * @tparam T is the type of elements.
     * @tparam ord is the storage order of the matrix, it can be Row_wise or Column_wise
     */

    template<typename T, StorageOrder ord>
    class Matrix{
        private:
            map<Index, T> data;
            vector<size_t> I, II; // I is the inner vector of indexes, II is the outer one in case the matrix is compressed
            vector<T> V; // V is the vector of the values of nnz elements in case the matrix is compressed
            size_t n_rows, n_cols, nnz_elem=0;
            bool compressed; 
            bool is_inbound(size_t i, size_t j) const;
            size_t calculate_nnz() const;

        public:
            Matrix(size_t r, size_t c): n_rows(r), n_cols(c), compressed(false){
                nnz_elem = calculate_nnz();
                resize(r,c);

            };
            void resize(size_t i, size_t j);
            void compress();
            void uncompress(); 
            bool is_compressed() const;
            T operator()(size_t i, size_t j) const; 
            T& operator()(size_t i, size_t j); 
            size_t get_n_rows() const;
            size_t get_n_cols() const;
            size_t get_nnz() const;
            void set_nnz(size_t nnz);
            T norm_one() const;
            T norm_infinity() const;
            T norm_frobenius() const;
            template<Norm n>
            T norm() const;
            void print_matrix() const;
            void print_compressed_matrix() const;
            void read_matrix(const string & filename);

            friend bool operator<(const Index & lhs, const Index & rhs);       

            template<typename U, StorageOrder order>
            friend vector<U> operator*(const Matrix<U,order> & mat, const vector<U> & vec);     

            template<typename U, StorageOrder order>
            friend Matrix<U,order> operator*(const Matrix<U,order> & matl, const Matrix<U,order> & matr);     

            
    };
}

#include "function_implementation.hpp"
#endif 