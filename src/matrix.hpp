//#ifndef matrix_hpp
//#define matrix_hpp


#include<iostream>
#include<vector>
#include<array>
#include<map>
#include <stdio.h>
#include <stdlib.h>

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
     * @brief Matrix a dynamic matrix class template.
     * @tparam T is the type of elements.
     * @tparam ord is the storage order of the matrix, it can be Row_wise or Column_wise
     */

    template<typename T, StorageOrder ord>
    class Matrix{
        private:
            using Index = array<size_t,2>;
            map<Index, T> data;
            vector<size_t> I, II; // I is the first vector of indexes, II is the second one in case the matrix is compressed
            vector<T> V; // V is the vector of the values of nnz elements in case the matrix is compressed
            size_t n_rows, n_cols, nnz_elem;
            bool compressed; 
            bool is_inbound(size_t i, size_t j) const;

        public:
            Matrix(size_t r, size_t c): n_rows(r), n_cols(c), compressed(false){};
            void resize(size_t i, size_t j);
            void compress();
            void uncompress(); 
            bool is_compressed() const;
            T operator()(size_t i, size_t j) const; 
            T operator()(size_t i, size_t j); 
            Matrix<T,ord> read_matrix(const string & filename) const;

            template<typename T, StorageOrder order>
            friend vector<T> operator*(const Matrix<T,order> & mat, const vector<T> & vec);      
            
    };

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
}


//#endif 