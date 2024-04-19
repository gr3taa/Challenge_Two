//#ifndef matrix_hpp
//#define matrix_hpp


#include<iostream>
#include<vector>
#include<array>
#include<map>
#include<Eigen/Dense>

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
            using Index = array<size_t,2>
            map<Index, T> data;
            vector<size_t> I, II; // I is the first vector of indexes, II is the second one in case the matrix is compressed
            vector<T> V; // V is the vector of the values of nnz elements in case the matrix is compressed
            size_t n_rows, n_cols, nnz_elem;
            bool compressed; 
            bool is_inbound(size_t i, size_t j) const;

        public:
            Matrix(size_t r, size_t c): n_rows(r), n_cols(c){};
            void resize(size_t i, size_t j);
            void compress();
            void uncompress(); 
            bool is_compressed() const;
            T operator()(size_t i, size_t j) const; 
            T operator()(size_t i, size_t j, const T & value); 
            
            template<typename U>
            friend vector<U> operator*(const vector<T> & vec); //const Matrix<U, Order>& mat

    };

}


//#endif 