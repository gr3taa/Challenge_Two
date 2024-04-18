#include"matrix.hpp";

namespace algebra{

/*!
 * @brief resizing the matrix.
 * @tparam T is the type of elements.
 * @tparam order is the storage order of the matrix.
 * @param i is the new size of rows.
 * @param j is the new size of colomns.
 */
template<typename T, StorageOrder order>
void Matrix<T,order>::resize(size_t i, size_t j){

}

/*!
 * @brief Converting the internal storage from a map to a compresses sparse row on comlumns depending on order.
 * @tparam T is the type of elements.
 * @tparam order is the storage order of the matrix.
 */
template<typename T, StorageOrder order>
void Matrix<T,order>::compress(){

}

/*!
 * @brief Bringing back the matrix to the uncompressed state.
 * @tparam T is the type of elements.
 * @tparam order is the storage order of the matrix.
 */
template<typename T, StorageOrder order>
void Matrix<T,order>::uncompress(){

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