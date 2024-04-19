#include"matrix.hpp";

namespace algebra{

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
    return i < n_rows && j < n_cols;
}

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
        if(!is_compressed){

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
    if(compre == true){
        cerr<<"The matrix is already compressed"<<endl;
        exit(1);
    }
    size_t a, b, c;
    if(order == "Row_wise"){
        a = 1;
        b = 0;
        c = n_rows;
    } else if(order == "Column_wise"){
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
    size_t x = -1;
    for (size_t i = 0; i < nnz_elem; ++i){
        V.push_back(data[i].second);
        II.push_back(data[i].first[a]);
        if(data[i].first[b]>x){
            II.push_back(i);
            ++x;
        }
    }
    I.push_back(c+R[0]);
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
    if(compre == false){
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
        if(k>=I[k]){
            ++t;
        }
        if(order == "Row_wise"){
            Index ij = {I[t],II[k]};
        } else if(order == "Row_wise"){
            Index ij = {II[k],I[t]};
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
 * @brief This function checks if the matrix is compressed or not.
 * @tparam T is the type of elements.
 * @tparam order is the storage order of the matrix.
 * @return true if the matrix is compressed.
 */
template<typename T, StorageOrder order>
bool Matrix<T,order>::is_compressed() const{
    return compressed;
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
    }
//TODO: compressed case
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
T Matrix<T,order>::operator()(size_t i, size_t j, const T & value){
    
}




}