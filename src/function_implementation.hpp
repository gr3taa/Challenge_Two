#include"matrix.hpp"

namespace algebra{

//01sparse matrix

    /**
    * @brief Check if i and j are whitin the range of the matrix.
    *
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

    /**
    * @brief This function calculates the number of non zero elements of the matrix.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @return The number of non zero elements of the matrix.
    */
    template<typename T, StorageOrder order>
    size_t Matrix<T,order>::calculate_nnz() const{
        return data.size();
    }

    /**
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
    /**
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

    /**
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
        size_t c;
        if(order == StorageOrder::Row_wise){
            c = n_rows;
        } else if(order == StorageOrder::Column_wise){
            c = n_cols;
        } else{
            cerr<<"Invalid Storage Order"<<endl;
            exit(1);
        }
        II.resize(nnz_elem);
        I.resize(c + 1);
        V.resize(nnz_elem);
        size_t x = 0, i = 0;
        I[0] =0;
        for (auto it = data.begin(); it != data.end(); ++it){
            
            V[i] = it->second;
            II[i] = it->first[1];
            if(it->first[0]>x){
                ++x;
                I[x] = i ;
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
        I[c] = nnz_elem;
        compressed = true;
        data.clear();
        
    }

        


    /**
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
        vector<size_t> R,C;
        size_t t=0;
        if constexpr (order == StorageOrder::Row_wise){
            for(size_t k = 0; k < nnz_elem; ++k){
                if(k == I[t+1]){
                    ++t;
                }
                Index ij = {t,II[k]};
                data[ij] = V[k];
            }
        }else if constexpr(order == StorageOrder::Column_wise){
            for(size_t k = 0; k < nnz_elem; ++k){
                if(k==I[t+1]){
                    ++t;
                }
                Index ij = {II[k],t};
                data[ij] = V[k];
            }
        } else{
            cerr<<"Invalid Storage Order"<<endl;
            exit(1);
        }
        I.clear();
        II.clear();
        V.clear();
        compressed = false;
        
    }

    /**
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
            Index ij;
            if(order == StorageOrder::Row_wise){
                ij = {i,j};
            } else if(order == StorageOrder::Column_wise){
                ij = {j,i};
            } else{
                cerr<<"Invalid Storage Order"<<endl;
                exit(1);
            }
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


    /**
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
            Index ij;
            if(order == StorageOrder::Row_wise){
                ij = {i,j};
            } else if(order == StorageOrder::Column_wise){
                ij = {j,i};
            } else{
                cerr<<"Invalid Storage Order"<<endl;
                exit(1);
            }
            ++nnz_elem;
            return data[ij];

        }else{
            cerr<<"The matrix is compressed"<<endl;
            exit(1);
        }
    }
    
    /**
    * @brief Get the number of rows of the matrix.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @return the number of rows of the matrix.
    */  
    template<typename T, StorageOrder order>
    size_t Matrix<T,order>::get_n_rows() const{
        return n_rows;
    }

    /**
    * @brief Get the number of columns of the matrix.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @return the number of columns of the matrix.
    */  
    template<typename T, StorageOrder order>
    size_t Matrix<T,order>::get_n_cols() const{
        return n_cols;
    }

    /**
    * @brief Get the number of non zero elements of the matrix.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @return the number of non zero elements of the matrix.
    */  
    template<typename T, StorageOrder order>
    size_t Matrix<T,order>::get_nnz() const{
        return nnz_elem;
    }

    /**
    * @brief Set the number of non zero elements of the matrix.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param nnz is the number of non zero elements.
    */ 
    template<typename T, StorageOrder order>
    void Matrix<T,order>::set_nnz(size_t nnz){
        nnz_elem = nnz;
    }

    /**
    * @brief this function computes the infinity norm.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param mat the matrix of which we want to calculate the norm.
    * @return the one norm.
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
        cout<<"The infinity norm is: "<<res<<endl;
        return res;
    }
    /**
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
        cout<<"The one norm is: "<<res<<endl;
        return res;
    }

    /**
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
    
    /**
    * @brief this function calculates the norm of the matrix.
    * @tparam U is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @tparam n is the type of norm.
    * @param mat is the matrix that we want to calculate the norm.
    * @return the norm of the matrix.
    */
    template<typename T, StorageOrder ord>
    template<Norm n>
    T Matrix<T,ord>::norm() const{
        if (n == Norm::One){
            return norm_one();
        } else if(n == Norm::Infinity){
            return norm_infinity();
        } else if(n == Norm::Frobenius){
            return norm_frobenius();
        } else{
            cerr<<"Error: Norm not recognized"<<endl;
            exit(1);
        }
    }

    template<typename T, StorageOrder order>
    void Matrix<T,order>::print_matrix() const{
        for (auto it = data.begin(); it != data.end(); ++it){
            cout<<"("<<it->first[0]<<", "<< it->first[1]<<") = " <<it->second<<" "<< endl;
        }
    }

    template<typename T, StorageOrder order>
    void Matrix<T,order>::print_compressed_matrix() const {
        if(order == StorageOrder::Row_wise){
            for(size_t i = 0; i < n_rows; i++){
                for(size_t j = I[i]; j < I[i+1]; j++){
                    cout<<"("<<i<<", "<<II[j]<<") = "<<V[j]<<endl;
                }
            }
        } else if(order == StorageOrder::Column_wise){
            for(size_t j = 0; j < n_cols; j++){
                for(size_t i = I[j]; i < I[j+1]; i++){
                    cout<<"("<<II[i]<<", "<<j<<") = "<<V[i]<<endl;
                }
            }
        } else{
            cerr<<"Error: Storage order not recognized"<<endl;
            exit(1);
        }
        for(size_t k = 0; k < nnz_elem; ++k){
            cout<<V[k]<<" - "<<II[k]<<endl;
            cout<<endl;
        }
        cout<<"I: "<<endl;
        for(size_t k = 0; k < n_rows+1; ++k){
            cout<<I[k]<<endl;
        }   
    }

    /**
    * @brief this function reads the matrix written in matrix market format.
    * @tparam T is the type of elements.
    * @tparam order is the storage order of the matrix.
    * @param filename is the name of the file containing the matrix.
    * @return the matrix.
    */
    template<typename T, StorageOrder order>
    void Matrix<T,order>::read_matrix(const string & filename) {
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
        iss>>n_rows>>n_cols>>nnz_elem;
        for(size_t k = 0; k < nnz_elem; ++k){
            size_t i,j;
            T val;
            file>>i>>j>>val;
            Index ij = {i,j};
            data[ij] = val;
        }
        set_nnz(nnz_elem);
    }

    
    /**
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

    /**
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

    /**
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

}