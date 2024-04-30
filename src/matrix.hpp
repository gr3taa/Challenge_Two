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
            vector<size_t> I, II; // I is the first vector of indexes, II is the second one in case the matrix is compressed
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
            T norm_one() const;
            T norm_infinity() const;
            T norm_frobenius() const;

            friend bool operator<(const Index & lhs, const Index & rhs);            
            
            template<typename U, StorageOrder order>
            friend Matrix<U,order> read_matrix(const string & filename); 

            template<typename U, StorageOrder order>
            friend vector<U> operator*(const Matrix<U,order> & mat, const vector<U> & vec);     

            template<typename U, StorageOrder order, Norm n>
            friend U norm(const Matrix<U,order> & mat);         
            
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
    * @brief this function reads the matrix written in matrix market format.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param filename is the name of the file containing the matrix.
    * @return the matrix.
    */
    template<typename T, StorageOrder order>
    Matrix<T,order> read_matrix(const string & filename) {
        ifstream file(filename);
        if(!file.is_open()){
            cerr<<"Error: file not found"<<endl;
            exit(1);
        }
        string line;
        getline(file,line);
        if(line != "%%MatrixMarket matrix coordinate real general"){
            cerr<<"Error: file not in Matrix Market format"<<endl;
            exit(1);
        }
        while(getline(file,line)){
            if(line[0] != '%'){
                break;
            }
        }
        istringstream iss(line);
        size_t n_rows, n_cols, nnz_elem;
        iss>>n_rows>>n_cols>>nnz_elem;
        Matrix<T,order> mat(n_rows,n_cols);
        
        for(size_t k = 0; k < nnz_elem; ++k){
            size_t i,j;
            T val;
            file>>i>>j>>val;
            mat(i,j) = val;
        }
        
        return mat;
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
    * @brief this function calculates the norm of the matrix.
    * @tparam U is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @tparam n is the type of norm.
    * @param mat is the matrix that we want to calculate the norm.
    * @return the norm of the matrix.
    */
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
    




}


#endif 