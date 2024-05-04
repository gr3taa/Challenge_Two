#include"matrix_class.hpp"

namespace algebra{

    
    template<typename T, StorageOrder order>
    bool Matrix<T,order>::is_inbound(size_t i, size_t j) const{
        return i < n_rows && j < n_cols;
    }

    
    template<typename T, StorageOrder order>
    size_t Matrix<T,order>::calculate_nnz() const{
        return data.size();
    }

    
    template<typename T, StorageOrder order>
    bool Matrix<T,order>::is_compressed() const{
        return compressed;
    }

    
    template<typename T, StorageOrder order>
    void Matrix<T,order>::resize(size_t i, size_t j){
        if(!is_compressed()){ 
            if(order == StorageOrder::Row_wise){
                n_rows = i;
                n_cols = j;
            } else if(order == StorageOrder::Column_wise){
                n_rows = j;
                n_cols = i;
            } else{
                cerr<<"Invalid Storage Order"<<endl;
                exit(1);
            }
        }
    }

    
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
        I[c] = nnz_elem;
        compressed = true;
        data.clear();
    }

    template<typename T, StorageOrder order>
    void Matrix<T,order>::uncompress(){
        if(compressed == false){
            cerr<<"The matrix is already uncompressed"<<endl;
            exit(1);
        }
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
            for(size_t j = 0; j < n_cols; ++j){
                for(size_t i = II[j]; i < II[j+1]; ++i){
                    data[{j,II[i]}] = V[i];
                }
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

    
    template<typename T, StorageOrder order>
    T Matrix<T,order>::operator()(size_t i, size_t j) const{
        if(!is_inbound(i,j)){
            cerr<<"The indexes are not in the range of the matrix"<<endl;
            exit(1);
        }
        if(!is_compressed()){
            Index ij = {0,0};
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
                return data.at(ij);
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
                for(size_t k = II[j]; k < II[j+1]; ++k){
                    if(I[k] == i){
                        return V[k];
                    }
                }
            } else{
                cerr<<"Invalid Storage Order"<<endl;
                exit(1);
            }
        }
    }

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
            cout<<"The matrix is compressed. You can only change the value of existing non-zeroelements."<<endl;
            if constexpr (order == StorageOrder::Row_wise) {
                for (std::size_t k = I[i]; k < I[i + 1]; ++k) {
                    if (II[k] == j) {
                        return V[k];
                    }
                }
            } else if constexpr (order == StorageOrder::Column_wise) {
                for (std::size_t k = II[j]; k < II[j + 1]; ++k) {
                    if (I[k] == i) {
                        return V[k];
                    }
                }
            }
        }
    }
    
    
    template<typename T, StorageOrder order>
    size_t Matrix<T,order>::get_n_rows() const{
        return n_rows;
    }

     
    template<typename T, StorageOrder order>
    size_t Matrix<T,order>::get_n_cols() const{
        return n_cols;
    }

    
    template<typename T, StorageOrder order>
    size_t Matrix<T,order>::get_nnz() const{
        return nnz_elem;
    }

    
    template<typename T, StorageOrder order>
    void Matrix<T,order>::set_nnz(size_t nnz){
        nnz_elem = nnz;
    }

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
            vector<T> sum(n_rows+1,0);
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
            vector<T> sum(n_cols+1,0);
                for(size_t k = 0; k < nnz_elem; ++k){
                    if(II[k] < n_cols +1)
                        sum[II[k]] += abs(V[k]);
                        else{
                            cerr<<"Error: The column index is out of range"<<endl;
                            exit(1);
                        }
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
    }

    
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
        if constexpr(order == StorageOrder::Row_wise){
            iss>>n_rows>>n_cols>>nnz_elem;
        } else if constexpr(order == StorageOrder::Column_wise){
            iss>>n_cols>>n_rows>>nnz_elem;
        } else{
            cerr<<"Error: Storage order not recognized"<<endl;
            exit(1);
        }
        for(size_t k = 0; k < nnz_elem; ++k){
            size_t i,j;
            T val;
            if constexpr(order == StorageOrder::Row_wise){
                file>>i>>j>>val;
            } else if constexpr(order == StorageOrder::Column_wise){
                file>>j>>i>>val;
            } else{
                cerr<<"Error: Storage order not recognized"<<endl;
                exit(1);
            }
            data[{i-1,j-1}] = val;
        }
        set_nnz(nnz_elem);
    }

    


}