#include"matrix_class.hpp"
namespace algebra{
        
        template<typename T, StorageOrder ord>
        bool operator<(const Index & lhs, const Index & rhs){
            if(ord == StorageOrder::Row_wise){
                return lhs[1] < rhs[1] || lhs[1] == rhs[1] && lhs[0] < rhs[0];
            }
            else if(ord == StorageOrder::Column_wise){
                return lhs[0] < rhs[0] || lhs[0] == rhs[0] && lhs[1] < rhs[1];
            } else{
                cerr<<"Error: Storage order not recognized"<<endl;
                exit(1);
            }
        }

        
        template<typename T, StorageOrder order>
        auto operator*(const Matrix<T,order> & mat, const vector<T> & vec){
            vector<common_type_t<T>> res(mat.n_rows);
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
                        for(size_t i = mat.I[j]; i < mat.I[j+1]; i++){
                            res[mat.II[i]] += mat.V[i]*vec[j];
                        }
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
            if(matl.is_compressed()){
                Matrix<U,order> res(matl.get_n_rows(), matr.get_n_cols());
                if(order == StorageOrder::Row_wise){
                    for(size_t i = 0; i < matl.get_n_rows(); i++){
                        for(size_t j = 0; j < matr.get_n_cols(); j++){
                            U sum = 0;
                            for(size_t k = 0; k < matr.get_n_rows(); k++){
                                for(size_t l = matl.I[i]; l < matl.I[i+1]; l++){
                                    if(matl.II[l] == k){
                                        for(size_t m = matr.I[k]; m < matr.I[k+1]; m++){
                                            if(matr.II[m] == j){
                                                sum += matl.V[l]*matr.V[m];
                                            }
                                        }
                                    }
                                }
                            }
                            res(i,j) = sum;
                        }
                    }
                } else if(order == StorageOrder::Column_wise){
                    for(size_t i = 0; i < matl.get_n_rows(); i++){
                        for(size_t j = 0; j < matr.get_n_cols(); j++){
                            U sum = 0;
                            for(size_t k = 0; k < matr.get_n_rows(); k++){
                                for(size_t l = matl.I[i]; l < matl.I[i+1]; l++){
                                    if(matl.II[l] == k){
                                        for(size_t m = matr.I[k]; m < matr.I[k+1]; m++){
                                            if(matr.II[m] == j){
                                                sum += matl.V[l]*matr.V[m];
                                            }
                                        }
                                    }
                                }
                            }
                            res(i,j) = sum;
                        }
                    }
                } else{
                    cerr<<"Error: Storage order not recognized"<<endl;
                    exit(1);
                }
                return res;
            } else if(!matl.is_compressed()){
                Matrix<U,order> res(matl.get_n_rows(), matr.get_n_cols());
                for(size_t i = 0; i < matl.get_n_rows(); ++i)
                    for(size_t j = 0; j < matr.get_n_cols(); ++j)
                        for(size_t k = 0; k < matr.get_n_rows(); ++k)
                            res(i,j) += matl(i,k) * matr(k,j);
                return res;
            } 
        }  
    }