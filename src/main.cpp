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
    Matrix<double, StorageOrder::Row_wise> mat(131,131);
    mat.read_matrix("lnsp_131.mtx");

    // matrix-vector product
    cout<<"-------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Matrix-vector product"<<endl;
    vector<double> vec(mat.get_n_cols(),2.0);
    Chrono timer_p;
    timer_p.start();
    vector<double> res = mat*vec;
    timer_p.stop();
    cout<<"Result of the matrix-vector product: "<<endl;
    print_vector(res);
    cout<<"Time for matrix-vector product: "<<timer_p.wallTime()<<" microseconds"<<endl;

    //compress
    cout<<"\n--------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Compression"<<endl;
    Chrono timer_c;
    timer_c.start();
    mat.compress();
    timer_c.stop();
    cout<<"Time for compression: "<<timer_c.wallTime()<<" microseconds"<<endl;
    
    //uncompress
    cout<<"\n--------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Uncompression"<<endl;    
    Chrono timer_u;
    timer_u.start();
    mat.uncompress();
    timer_u.stop();
    cout<<"Time for uncompression: "<<timer_u.wallTime()<<" microseconds"<<endl;

    

    //norm
    cout<<"\n--------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Norm"<<endl;
    Chrono timer_n;
    timer_n.start();
    double norm = mat.norm(Norm::Frobenius);
    timer_n.stop();
    cout<<"Time for norm calculation: "<<timer_n.wallTime()<<" microseconds"<<endl;


    return 0;
};

// lineareAlgebra