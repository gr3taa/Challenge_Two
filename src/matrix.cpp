#include"matrix.hpp"

namespace algebra{

    template class Matrix<double, algebra::StorageOrder::Row_wise>;
    template class Matrix<double, algebra::StorageOrder::Column_wise>;

    /*!
    * @brief Check if i and j are whitin the range of the matrix.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param i is the index of the row.
    * @param j is the index of the column.
    * @return true if i and j are whitin the range of the matrix, false otherwise.
    */
    template<typename T, StorageOrder order>
    bool Matrix<T,order>::is_inbound(size_t i, size_t j) const{
        return i <= n_rows && j <= n_cols;
    }

    /*!
    * @brief This function calculates the number of non zero elements of the matrix.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @return The number of non zero elements of the matrix.
    */
    template<typename T, StorageOrder order>
    size_t Matrix<T,order>::calculate_nnz() const{
        return data.size();
    }

    /*!
    * @brief This function checks if the matrix is compressed or not.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @return true if the matrix is compressed.
    */
    template<typename T, StorageOrder order>
    bool Matrix<T,order>::is_compressed() const{
        return compressed;
    }

    //TOREview
    /*!
    * @brief resizing the matrix.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param i is the new size of rows.
    * @param j is the new size of colomns.
    */
    template<typename T, StorageOrder order>
    void Matrix<T,order>::resize(size_t i, size_t j){
        if(i < n_rows || j < n_cols){
            if(!is_compressed()){
                for(auto it = data.begin(); it != data.end();){
                    if(it->first[0] >= i || it->first[1] >= j){
                        it = data.erase(it);
                    } else{
                        ++it;
                    }
                }
            } else{
                for(size_t k = 0; k < nnz_elem; ++k){
                    if(II[k] >= i || I[k] >= j){
                        II.erase(II.begin() + k);
                        V.erase(V.begin() + k);
                        --k;
                        --nnz_elem;
                    }
                }
                I.clear();
                I.push_back(0);
                size_t x = -1;
                for(size_t k = 0; k < nnz_elem; ++k){
                    if(II[k] > x){
                        I.push_back(k);
                        ++x;
                    }
                }
            }
        }
        //TODO: eliminare gli elementi out of range
        n_rows = i;
        n_cols = j;
    }

    /*!
    * @brief Converting the internal storage from a map to a compresses sparse row on comlumns depending on order.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    */
    template<typename T, StorageOrder order>
    void Matrix<T,order>::compress(){ 
        if(compressed == true){
            cerr<<"The matrix is already compressed"<<endl;
            exit(1);
        }
        size_t a, b, c;
        if(order == StorageOrder::Row_wise){
            a = 1;
            b = 0;
            c = n_rows;
        } else if(order == StorageOrder::Column_wise){
            a = 0;
            b = 1;
            c = n_cols;
        } else{
            cerr<<"Invalid Storage Order"<<endl;
            exit(1);
        }
        II.resize(nnz_elem);
        I.resize(c + 1);
        V.resize(nnz_elem);
        size_t x = -1, i = 0;
        for (auto it = data.begin(); it != data.end(); ++it){
            V.push_back(it->second);
            II.push_back(it->first[a]);
            if(it->first[b]>x){
                II.push_back(i);
                ++x;
            }
            ++i;
        }
/*
        for (size_t i = 0; i < nnz_elem; ++i){
            V.push_back(data[i].second);
            II.push_back(data[i].first[a]);
            if(data[i].first[b]>x){
                II.push_back(i);
                ++x;
            }
        }*/
        I.push_back(c+I[0]);
        compressed = true;
        data.clear();
        
    }

        


    /*!
    * @brief Bringing back the matrix to the uncompressed state.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    */
    template<typename T, StorageOrder order>
    void Matrix<T,order>::uncompress(){
        if(compressed == false){
            cerr<<"The matrix is already uncompressed"<<endl;
            exit(1);
        }
        /*vector<size_t> R,C;
        size_t t=0;
        if(order == "Row_wise"){
            for(size_t k = 0; k < nnz_elem; ++k){
            if(k>=I[k]){
                ++t;
            }
            Index ij = {I[t],II[k]};
            data[ij] = V[k];
            }
        }else if(order == "Column_wise"){
            for(size_t k = 0; k < nnz_elem; ++k){
            if(k>=I[k]){
                ++t;
            }
            Index ij = {II[k],I[t]};
            data[ij] = V[k];
        }
        } else{
            cerr<<"Invalid Storage Order"<<endl;
            exit(1);
        }*/
        size_t t = 0;
        for(size_t k = 0; k < nnz_elem; ++k){
            Index ij;
            if(k>=I[k]){
                ++t;
            }
            if(order == StorageOrder::Row_wise){
                ij = {I[t],II[k]};
            } else if(order == StorageOrder::Column_wise){
                ij = {II[k],I[t]};
            } else{
            cerr<<"Invalid Storage Order"<<endl;
            exit(1);
            }
            data[ij] = V[k];
            }
        I.clear();
        II.clear();
        V.clear();
        compressed = false;
        
    }

    /*!
    * @brief Overloading call operator
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param i is the index of the row.
    * @param j is the index of the column.
    * @return 0 if the element is within the range of the matrix but is not present.
    */
    template<typename T, StorageOrder order>
    T Matrix<T,order>::operator()(size_t i, size_t j) const{
        if(!is_inbound(i,j)){
            cerr<<"The indexes are not in the range of the matrix"<<endl;
            exit(1);
        }
        if(!is_compressed()){
            Index ij = {i,j};
            auto it = data.find(ij);
            if(it!= data.end())
                return it->second;
            else
                return 0;
        }else{
            if(order == StorageOrder::Row_wise){
                for(size_t k = I[i]; k < I[i+1]; ++k){
                    if(II[k] == j){
                        return V[k];
                    }
                }
            } else if(order == StorageOrder::Column_wise){
                for(size_t k = I[j]; k < I[j+1]; ++k){
                    if(II[k] == i){
                        return V[k];
                    }
                }
            } else{
                cerr<<"Invalid Storage Order"<<endl;
                exit(1);
            }
        }
    }


    /*!
    * @brief Overloading call operator.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param i is the index of the row.
    * @param j is the index of the column.
    * @return 
    */
    template<typename T, StorageOrder order>
    T& Matrix<T,order>::operator()(size_t i, size_t j){
        if(!is_inbound(i,j)){
            cerr<<"The indexes are not in the range of the matrix"<<endl;
            exit(1);
        }
        if(!is_compressed()){
            Index ij = {i,j};
            return data[ij];
        }else{
            cerr<<"The matrix is compressed"<<endl;
            exit(1);
        }
    }
    /*!
    * @brief Get the number of rows of the matrix.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @return the number of rows of the matrix.
    */  
    template<typename T, StorageOrder order>
    size_t Matrix<T,order>::get_n_rows() const{
        return n_rows;
    }

    /*!
    * @brief Get the number of columns of the matrix.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @return the number of columns of the matrix.
    */  
    template<typename T, StorageOrder order>
    size_t Matrix<T,order>::get_n_cols() const{
        return n_cols;
    }

    /*!
    * @brief Get the number of non zero elements of the matrix.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @return the number of non zero elements of the matrix.
    */  
    template<typename T, StorageOrder order>
    size_t Matrix<T,order>::get_nnz() const{
        return nnz_elem;
    }

    /*!
    * @brief this function computes the one norm.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param mat the matrix of which we want to calculate the norm.
    * @return the one norm.
    */
    template <typename T, StorageOrder order>
    T Matrix<T,order>::norm_one() const{
        T res = 0;
        if(!is_compressed()){
            for(size_t j = 0; j < n_cols; ++j){
                T sum = 0;
                for(size_t i = 0; i < n_rows; ++i){
                    sum += abs((*this)(i,j));
                }
                if(sum > res)
                    res = sum;
            }
        }else if (order == StorageOrder::Column_wise){
            vector<T> sum(n_rows,0);
            for(size_t k = 0; k < nnz_elem; ++k){
                sum[II[k]] += abs(V[k]);
            }
            res = *max_element(sum.begin(),sum.end());
        }else if (order == StorageOrder::Row_wise){
            for(size_t i = 0; i < n_rows; ++i){
                T sum = 0;
                for(size_t k = I[i]; k < I[i+1]; ++k){
                    sum += abs(V[k]);
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
    T Matrix<T,order>::norm_infinity() const{
        T res = 0;
        if(!is_compressed()){
            for(size_t i = 0; i < n_rows; ++i){
                T sum = 0;
                for(size_t j = 0; j < n_cols; ++j){
                    sum += abs((*this)(i,j));
                }
                if(sum > res)
                    res = sum;
            }
            
        }else if (order == StorageOrder::Row_wise){
            vector<T> sum(n_rows,0);
                for(size_t k = 0; k < nnz_elem; ++k){
                    sum[II[k]] += abs(V[k]);
                }
                res = *max_element(sum.begin(),sum.end());
        }else if (order == StorageOrder::Column_wise){
                for(size_t j = 0; j < n_cols; ++j){
                    T sum = 0;
                    for(size_t k = I[j]; k < I[j+1]; ++k){
                        sum += abs(V[k]);
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
    T Matrix<T,order>::norm_frobenius() const{
        T res = 0;
        if(!is_compressed()){
            for(auto it = data.begin(); it != data.end(); ++it){
                res += pow(it->second,2);
                
            }
        }else{
            for(size_t k = 0; k < nnz_elem; ++k){
                res += pow(V[k],2);
            }
        }
        res = sqrt(res);
        cout<<"The frobenius norm is: "<<res<<endl;
        return res;
    }
    

}