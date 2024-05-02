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
            //template<Norm n>
            T norm(Norm n) const;
            void print_matrix() const;
            void print_compressed_matrix() const;
            void read_matrix(const string & filename);

            friend bool operator<(const Index & lhs, const Index & rhs);       

            template<typename U, StorageOrder order>
            friend vector<U> operator*(const Matrix<U,order> & mat, const vector<U> & vec);     

            template<typename U, StorageOrder order>
            friend Matrix<U,order> operator*(const Matrix<U,order> & matl, const Matrix<U,order> & matr);     

/*
            template<typename U, StorageOrder order, Norm n>
            friend U norm(const Matrix<U,order> & mat);       */  
            
    };

    /*!
    * @brief Comparison operator according to the storage order.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param lhs is the first index.
    * @param rhs is the second index.
    * @return true if lhs is less than rhs, false otherwise.
    */
    template<typename T, StorageOrder ord>
    bool operator<(const Index & lhs, const Index & rhs){
        if(ord == StorageOrder::Column_wise){
            return lhs[1] < rhs[1] || (lhs[1] == rhs[1] && lhs[0] < rhs[0]);
        }
        else if(ord == StorageOrder::Row_wise){
            return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]);
        } else{
            cerr<<"Error: Storage order not recognized"<<endl;
            exit(1);
        }
    }

    /*!
    * @brief this function does the matrix-vector product.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param mat is the matrix that we want to multiply.
    * @param vec is the vector that we want to multiply.
    * @return the result of the matrix-vector product.
    */
    template<typename T, StorageOrder order>
    vector<T> operator*(const Matrix<T,order> & mat, const vector<T> & vec){
        vector<T> res(mat.n_rows);
        if(mat.is_compressed()){
            if(order == StorageOrder::Row_wise){
                for(size_t i = 0; i < mat.n_rows; i++){
                    T sum = 0;
                    for(size_t j = mat.I[i]; j < mat.I[i+1]; j++){
                        sum += mat.V[j]*vec[mat.II[j]];
                    }
                    res[i] = sum;
                }
            }else if(order == StorageOrder::Column_wise){
                for(size_t j = 0; j < mat.n_cols; j++){
                    T sum = 0;
                    for(size_t i = mat.I[j]; i < mat.I[j+1]; i++){
                        sum += mat.V[i]*vec[mat.II[i]];
                    }
                    res[j] = sum;
                }
            } else{
                cerr<<"Error: Storage order not recognized"<<endl;
                exit(1);
            }
        } else{
            for(size_t i = 0; i < mat.n_rows; i++){
                T sum = 0;
                for(size_t j = 0; j < mat.n_cols; j++){
                    sum += mat(i,j)*vec[j];
                }
                res[i] = sum;
            }
        }
        return res;
    }

    /*!
    * @brief this function does the matrix-matrix  product.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param mat is the matrix that we want to multiply.
    * @param mat1 is the matrix that we want to multiply.
    * @return the result of the matrix-matrix product.
    */
    template<typename U, StorageOrder order>
    Matrix<U,order> operator*(const Matrix<U,order> & matl, const Matrix<U,order> & matr){
        if(matl.get_n_cols() != matr.get_n_rows()){
            cerr<<"Error: The number of columns of the left matrix must be equal to the number of rows of the right matrix"<<endl;
            exit(1);
        }
        if(matl.is_compressed() != matr.is_compressed()){
            cerr<<"Error: The matrices must have the same storage order"<<endl;
            exit(1);
        }
        if(!matl.is_compressed()){
            Matrix<U,order> res(matl.get_n_rows(), matr.get_n_cols());
            for(size_t i = 0; i < matl.get_n_rows(); ++i)
                for(size_t j = 0; j < matr.get_n_cols(); ++j)
                    for(size_t k = 0; k < matr.get_n_rows(); ++k)
                        res(i,j) += matl(i,k) * matr(k,j);
            return res;
        } /*else if (order == StorageOrder::Row_wise){
                    Matrix<U,order> res(matl.get_n_rows(),matr.get_n_cols());
                    vector<U> vec(matr.get_n_rows());
                    for(size_t j = 0; j < matr.get_n_cols(); ++j){
                        // jth column of matr
                        vector<U> vec = get_col(matr,j);
                        for(size_t i = 0; i < matl.get_n_rows(); ++i){
                            size_t row_starts = matl.I[i];
                            res(i,j) = inner_product(vec.begin(), vec.end(), matl.V.cbegin() + row_starts, U(0));
                        }   
                    }
                    return res;
        }*/
    }  
/*
    template<typename T, StorageOrder ord>
    vector<T> get_col(const Matrix<T,ord> & mat, size_t j){
        if(order == StorageOrder::Column_wise || !is.compress()){
            cerr<<"Error: The matrix must be row-wise compressed"<<endl;
            exit(1);
        }
        vector<T> res(mat.get_n_rows());
        for(size_t i = 0; i < mat.get_n_rows(); i++){
            size_t row_start = mat.I[i];
            size_t row_end = mat.I[i+1];
            auto it = lower_bound(mat.II.begin() + row_start, mat.II.begin() + row_end, j);
            if(it != mat.II.begin() + row_end && *it == j){
                res[i] = mat.V[it - mat.II.begin()];
            } else{
                res[i] = 0;
            }
        }
        return res;
    }*/


    /*!
    * @brief this function calculates the norm of the matrix.
    * @tparam U is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @tparam n is the type of norm.
    * @param mat is the matrix that we want to calculate the norm.
    * @return the norm of the matrix.
    *//*
    template<typename U, StorageOrder order, Norm n>
    U norm(const Matrix<U,order> & mat){
        if(n == Norm::One){
            return mat.norm_one();
        } else if(n == Norm::Infinity){
            return mat.norm_infinity();
        } else if(n == Norm::Frobenius){
            return mat.norm_frobenius();
        } else{
            cerr<<"Error: Norm not recognized"<<endl;
            exit(1);
        }
    }
    
*/



}


#endif 