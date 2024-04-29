#include"matrix.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

namespace algebra{
    
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

    /**
     * @brief Enum defining the type of norm.
     */
    enum class Norm{
        One,
        Infinity,
        Frobenius
    };

    /*!
    * @brief this function computes the one norm.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param mat the matrix of which we want to calculate the norm.
    * @return the one norm.
    */
    template <typename T, StorageOrder order>
    T norm_one(const Matrix<T,order> & mat){
        double res = 0;
        if(!mat.is_compressed()){
            for(size_t j = 0; j < mat.get_n_cols(); ++j){
                T sum = 0;
                for(size_t i = 0; i < mat.get_n_rows(); ++i){
                    sum += abs(mat(i,j));
                }
                if(sum > res)
                    res = sum;
            }
        } else if (order == StorageOrder::Column_wise){
                vector<T> sum(mat.get_n_rows(),0);
                for(size_t k = 0; k < mat.get_nnz(); ++k){
                    sum(mat.II[k]) += abs(mat.V[k]);
                }
                res = *max_element(sum.begin(),sum.end());
        } else if (order == StorageOrder::Row_wise){
                for(size_t i = 0; i < mat.get_n_rows(); ++i){
                    T sum = 0;
                    for(size_t k = mat.I[i]; k < mat.I[i+1]; ++k){
                        sum += abs(mat.V[k]);
                    }
                    if(sum > res)
                        res = sum;  
                }
        } else{
            cerr<<"Error: Storage order not recognized"<<endl;
            exit(1);
        }
        cout<<"The one norm is: "<<res<<endl;
        return res;
    }
    
    /*!
    * @brief this function computes the infinity norm.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param mat the matrix of which we want to calculate the norm.
    * @return the infinity norm.
    */
    template <typename T, StorageOrder order>
    T norm_infinity(const Matrix<T,order> & mat){
        T res = 0;
        if(!mat.is_compressed()){
            for(size_t i = 0; i < mat.get_n_rows(); ++i){
                T sum = 0;
                for(size_t j = 0; j < mat.get_n_cols(); ++j){
                    sum += abs(mat(i,j));
                }
                if(sum > res)
                    res = sum;
            }
            
        }else if (order == StorageOrder::Row_wise){
            vector<T> sum(mat.get_n_rows(),0);
                for(size_t k = 0; k < mat.get_nnz(); ++k){
                    sum(mat.II[k]) += abs(mat.V[k]);
                }
                res = *max_element(sum.begin(),sum.end());
        }else if (order == StorageOrder::Column_wise){
                for(size_t j = 0; j < mat.get_n_cols(); ++j){
                    T sum = 0;
                    for(size_t k = mat.I[j]; k < mat.I[j+1]; ++k){
                        sum += abs(mat.V[k]);
                    }
                    if(sum > res)
                        res = sum;  
                }
        } else{
            cerr<<"Error: Storage order not recognized"<<endl;
            exit(1);
        }
        cout<<"The infinity norm is: "<<res<<endl;
        return res;
    }

    /*!
    * @brief this function computes the frobenius norm.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param mat the matrix of which we want to calculate the norm.
    * @return the frobenius norm.
    */
    template <typename T, StorageOrder order>
    T norm_frobenius(const Matrix<T,order> & mat){
        T res = 0;
        if(!mat.is_compressed()){
            for(auto it = mat.data.begin(); it != mat.data.end(); ++it){
                res += pow(it->second,2);
            }
        }
        else{
            T res = 0;
            for(size_t k = 0; k < mat.get_nnz(); ++k){
                res += pow(mat.V[k],2);
            }
        }
        cout<<"The Frobenius norm is: "<<sqrt(res)<<endl;
        return sqrt(res);
    }

    /*!
    * @brief this function computes the one or the infinity or the frobenius norm according to the template value.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @tparam n is the type of norm.
    * @param mat the matrix of which we want to calculate the norm.
    * @return the norm.
    */
    template<typename T, StorageOrder order, Norm n>
    T norm(const Matrix<T,order> & mat){
        T res = 0;
        if(n == Norm::One){
            res = norm_one<T,order>(mat);
        } else if(n == Norm::Infinity){
            res = norm_infinity<T,order>(mat);
        } else if(n == Norm::Frobenius){
            res = norm_frobenius<T,order>(mat);
        } else{
            cerr<<"Error: Norm not recognized"<<endl;
            exit(1);
        }
        return res;
    }
}
