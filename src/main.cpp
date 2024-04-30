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
    Matrix<double, StorageOrder::Column_wise> mat = read_matrix<double,StorageOrder::Column_wise>("lnsp_131.mtx");
    
    // matrix-vector product
    cout<<"Matrix-vector product---------------------------------------------------------------------------------------------"<<endl;
    vector<double> vec(mat.get_n_cols(),2.0);
    Chrono timer_p;
    timer_p.start();
    vector<double> res = mat*vec;
    timer_p.stop();
    cout<<"Time for matrix-vector product: "<<timer_p.wallTime()<<" microseconds"<<endl;
    cout<<"Result of the matrix-vector product: "<<endl;
    print_vector(res);

    //compress
    cout<<"Compression------------------------------------------------------------------------------------------------------"<<endl;
    Chrono timer_c;
    timer_c.start();
    mat.compress();
    timer_c.stop();
    cout<<"Time for compression: "<<timer_c.wallTime()<<" microseconds"<<endl;

    //uncompress
    cout<<"Uncompression----------------------------------------------------------------------------------------------------"<<endl;
    Chrono timer_u;
    timer_u.start();
    mat.uncompress();
    timer_u.stop();
    cout<<"Time for uncompression: "<<timer_u.wallTime()<<" microseconds"<<endl;


    return 0;
};

// lineareAlgebra