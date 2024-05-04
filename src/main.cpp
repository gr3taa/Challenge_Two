#include "matrix_class.hpp"
#include "chrono.hpp"
#include <iostream>
using namespace std;
using namespace algebra;
using namespace Timings;

void print_vector(const vector<double>& vec){
    for(auto i: vec){
        cout<<i<<"\n";
    }
    cout<<endl;
}


int main(){
    ifstream file("file/input.txt");
    string order, norm_type;
    if(!file.is_open()){
        cout<<"File not found"<<endl;
        return 0;
    }
    getline(file, order);
    getline(file, norm_type);
    Matrix<double, StorageOrder::Row_wise> mat(0,0);
    if(order == "Column_wise"){
        Matrix<double, StorageOrder::Column_wise> mat(0,0);
    } else if(order != "Row_wise" ){
        cout<<"Invalid storage order"<<endl;
        return 0;
    } 
    
    mat.read_matrix("file/lnsp_131.mtx");
    vector<double> vec(mat.get_n_cols(),2.0);

    // matrix-vector product (uncompressed version)
    cout<<"-------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Matrix-vector product (uncompressed version)"<<endl;
    Chrono timer_pu;
    timer_pu.start();
    vector<double> resp = algebra::operator*(mat, vec);
    timer_pu.stop();
    cout<<"Result of the matrix-vector product: "<<endl;
    print_vector(resp);
    cout<<"Time for matrix-vector product(uncompressed version): "<<timer_pu.wallTime()<<" microseconds"<<endl;

    //norm (uncompressed version)
    cout<<"\n--------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Norm"<<endl;
    Chrono timer_nu;
    timer_nu.start();
    double norm_u;
    if(norm_type == "One"){
        norm_u = mat.norm<Norm::One>();
    } else if(norm_type == "Infinity"){
        norm_u = mat.norm<Norm::Infinity>();
    } else if(norm_type == "Frobenius"){
        norm_u = mat.norm<Norm::Frobenius>();
    } else {
        cout<<"Invalid norm type"<<endl;
        return 0;
    }
    timer_nu.stop();
    cout<<"Time for norm calculation (uncompressed version): "<<timer_nu.wallTime()<<" microseconds"<<endl;

    //compress
    cout<<"\n--------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Compression"<<endl;
    Chrono timer_c;
    timer_c.start();
    mat.compress();
    timer_c.stop();
    cout<<"Time for compression: "<<timer_c.wallTime()<<" microseconds"<<endl;
    
    // matrix-vector product (compressed version)
    cout<<"-------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Matrix-vector product (compressed version)"<<endl;
    Chrono timer_p;
    timer_p.start();
    vector<double> res = algebra::operator*(mat, vec);
    timer_p.stop();
    cout<<"Result of the matrix-vector product: "<<endl;
    print_vector(res);
    cout<<"Time for matrix-vector product(compressed version): "<<timer_p.wallTime()<<" microseconds"<<endl;

    //norm (compressed version)
    cout<<"\n--------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Norm (compressed version)"<<endl;
    //mat.print_compressed_matrix();
    Chrono timer_n;
    timer_n.start();
    double norm ;
    if (norm_type == "One"){
        norm = mat.norm<Norm::One>();
    } else if(norm_type == "Infinity"){
        norm = mat.norm<Norm::Infinity>();
    } else if(norm_type == "Frobenius"){
        norm = mat.norm<Norm::Frobenius>();
    } else {
        cout<<"Invalid norm type"<<endl;
        return 0;
    }
    timer_n.stop();
    cout<<"Time for norm calculation (compressed version): "<<timer_n.wallTime()<<" microseconds"<<endl;


    //uncompress
    cout<<"\n--------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Uncompression"<<endl;    
    Chrono timer_u;
    timer_u.start();
    mat.uncompress();
    timer_u.stop();
    cout<<"Time for uncompression: "<<timer_u.wallTime()<<" microseconds"<<endl;


    
    return 0;
};

// lineareAlgebra