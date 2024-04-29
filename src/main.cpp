#include "helper.hpp"
#include "chrono.hpp"
#include <iostream>
using namespace std;
using namespace algebra;
using namespace Timings;

int main(){
    Matrix<double, StorageOrder::Row_wise> mat = read_matrix<double, StorageOrder::Row_wise>("Insp_131.mtx");

    vector<double> vec(mat.get_n_cols(),2.0);
    Chrono timer;
    timer.start();
    vector<double> res = mat*vec;
    timer.stop();
    cout<<"Time for matrix-vector product: "<<timer.wallTime()<<" microseconds"<<endl;


    return 0;
};

// lineareAlgebra