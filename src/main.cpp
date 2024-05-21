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

    // test with small matrices
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    std::cout << "small matrices - row" << std::endl;
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    
    Matrix<double, StorageOrder::Row_wise> M(3,3);
    M(0,0) = 2.0;
    M(1,1) = 3.0;
    M(2,2) = 4.0;

    std::vector<double> V = {1.0,1.0,1.0};
    
    M.norm<Norm::One>();
    M.norm<Norm::Infinity>();
    M.norm<Norm::Frobenius>();

    std::cout << "product: ";
    print_vector(M*V);

    M.print_matrix();
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    M.compress();
    M.print_compressed_matrix();
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    M.uncompress();
    M.print_matrix();
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    M.resize(2,2);
    M.print_matrix();
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    M.resize(5,5);
    M.print_matrix();
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
 
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    std::cout << "small matrices - col" << std::endl;
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    
    Matrix<double, StorageOrder::Column_wise> Mc(3,3);
    Mc(0,0) = 2.0;
    Mc(1,1) = 3.0;
    Mc(2,2) = 4.0;

    Mc.norm<Norm::One>();
    Mc.norm<Norm::Infinity>();
    Mc.norm<Norm::Frobenius>();

    std::cout << "product: ";
    print_vector(Mc*V);

    Mc.print_matrix();
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    Mc.compress();
    Mc.print_compressed_matrix();
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    Mc.uncompress();
    Mc.print_matrix();
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    Mc.resize(2,2);
    Mc.print_matrix();
    std::cout << "-----------------------------------------------------------------------------------" << std::endl; 
    Mc.resize(5,5);
    Mc.print_matrix();
    std::cout << "-----------------------------------------------------------------------------------" << std::endl; 

    // end test with small matrices

    ifstream file("file/input.txt");
    string order, norm_type;
    if(!file.is_open()){
        cout<<"File not found"<<endl;
        return 0;
    }
    getline(file, order);
    getline(file, norm_type);

    // Here you don't need to store both types of matrices
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
    //double norm_u;       // no need for storing the output in a variable
    if(norm_type == "One"){
        mat.norm<Norm::One>();
    } else if(norm_type == "Infinity"){
        mat.norm<Norm::Infinity>();
    } else if(norm_type == "Frobenius"){
        mat.norm<Norm::Frobenius>();
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
    //double norm ;       // no need for storing the output in a variable
    if (norm_type == "One"){
        mat.norm<Norm::One>();
    } else if(norm_type == "Infinity"){
        mat.norm<Norm::Infinity>();
    } else if(norm_type == "Frobenius"){
        mat.norm<Norm::Frobenius>();
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
