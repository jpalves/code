#include <boost/python/numpy.hpp>
#include <boost/python.hpp>
//#include <numpy/arrayobject.h>
#include <iostream>
#include "mat.h"
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>

namespace p  = boost::python;
namespace np = boost::python::numpy;

void printM(matriz <double> m){
    std::cout.precision(3);
    
    for (int i = 0; i < m.linhas(); ++i) {
        for (int j = 0; j < m.colunas(); ++j) {
            		std::cout << m[i][j] << "\t";
        	}
        	std::cout <<std::endl;
    	}std::cout << std::endl;
}

/// @brief 
/// @param A 
/// @param  
/// @return 
matriz <double> pinv(p::object r){
    np::ndarray uu = p::extract<np::ndarray>(r[0]);
    np::ndarray s  = p::extract<np::ndarray>(r[1]);
    np::ndarray vv = p::extract<np::ndarray>(r[2]);
    matriz <double> S(s.shape(0),s.shape(0),0.0);
    matriz <double> U(uu.shape(0),uu.shape(1));
    matriz <double> V(vv.shape(0),vv.shape(1));
    
    for(int i = 0; i < uu.shape(0); i++)
        for(int j = 0; j < uu.shape(1); j++)
            U[i][j] = reinterpret_cast<double*>(uu.get_data())[i+j*uu.shape(0)];
    
    for(int i = 0; i < vv.shape(0); i++)
        for(int j = 0; j < vv.shape(1); j++)
            V[i][j] = reinterpret_cast<double*>(vv.get_data())[i+j*vv.shape(0)];

    for(int i = 0; i < S.linhas(); i++)
        S[i][i] =  1/reinterpret_cast<double*>(s.get_data())[i];
    
    std::cout << S << std::endl;
    std::cout << V << std::endl;
    std::cout << U << std::endl;
    return V.transposta()*S*U.transposta();
}   

int main(int argc, char *argv[]){
    Py_Initialize();
    np::initialize();
    p::object scipy = p::import("scipy.linalg");
    p::object svd = scipy.attr("svd");
    double data[] = {0, 2, 3, 
                     4, 0, 0,
                     7, 8 ,9};

    np::ndarray A = np::from_data(data, np::dtype::get_builtin<double>(),
                                  p::make_tuple(9), p::make_tuple(sizeof(double)), p::object());
  
    matriz<double> INV = pinv(svd(A.reshape(p::make_tuple(3, 3))));
    //std::cout << std::setprecision(3) << std::fixed;
    //printM(INV);
    std::cout << INV << std::endl;
    std::cout << INV.transposta() << std::endl;

    return 0;
}