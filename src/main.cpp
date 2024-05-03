#include "matrix.hpp"
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
    Matrix<double, StorageOrder::Column_wise> mat(0,0);
    mat.read_matrix("lnsp_131.mtx");
    vector<double> vec(mat.get_n_cols(),2.0);

    // matrix-vector product (uncompressed version)
    cout<<"-------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Matrix-vector product (uncompressed version)"<<endl;
    Chrono timer_pu;
    timer_pu.start();
    vector<double> resp = mat*vec;
    timer_pu.stop();
    cout<<"Result of the matrix-vector product: "<<endl;
    print_vector(resp);
    cout<<"Time for matrix-vector product(uncompressed version): "<<timer_pu.wallTime()<<" microseconds"<<endl;

    //norm (uncompressed version)
    cout<<"\n--------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Norm"<<endl;
    Chrono timer_nu;
    timer_nu.start();
    double norm_u = mat.norm<Norm::One>();
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
    vector<double> res = mat*vec;
    timer_p.stop();
    cout<<"Result of the matrix-vector product: "<<endl;
    print_vector(res);
    cout<<"Time for matrix-vector product(compressed version): "<<timer_p.wallTime()<<" microseconds"<<endl;

    //norm (compressed version)
    cout<<"\n--------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Norm"<<endl;
    Chrono timer_n;
    timer_n.start();
    double norm = mat.norm<Norm::One>();
    timer_n.stop();
    cout<<norm_u<<"Time for norm calculation (compressed version): "<<timer_n.wallTime()<<" microseconds"<<endl;


    //uncompress
    cout<<"\n--------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Uncompression"<<endl;    
    Chrono timer_u;
    timer_u.start();
    mat.uncompress();
    timer_u.stop();
    cout<<"Time for uncompression: "<<timer_u.wallTime()<<" microseconds"<<endl;


    
/*
    Matrix<int,StorageOrder::Row_wise> A(3,3);
    Matrix<int, StorageOrder::Row_wise> B(3,2);
    A(0,0)=1;
    //A(0,1)=2;
    A(0,2)=2;
    A(1,0)=3;
    A(1,1)=4;
    //A(1,2)=6;
    //A(2,0)=7;
    //A(2,1)=8;
    A(2,2)=9;
    //A.print_matrix();
    cout<<"norm un"<<A.norm<Norm::Infinity>()<<endl;
    A.compress();
    cout<<"norm c"<<A.norm<Norm::Infinity>()<<endl;
*/
    /*B(0,0)=2;
    B(1,1)=5;
    B(2,0)=6;
    A.compress();
    B.compress();
    cout<<"Matrix A"<<endl;
    A.print_compressed_matrix();
    cout<<"Matrix B"<<endl;
    B.print_compressed_matrix();
    Matrix<int,StorageOrder::Row_wise> r = A * B;
    cout<<"Matrix r"<<endl;
    r.print_matrix();
*/
    
    return 0;
};

// lineareAlgebra