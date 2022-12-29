#include <iostream>
#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>

#include "common.hpp"
#include "mat.h"

template<typename T> void printM(matriz <T> m){
    std::cout.precision(3);
    
    for (int i = 0; i < m.linhas(); ++i) {
        for (int j = 0; j < m.colunas(); ++j) {
            		std::cout << m[i][j] << "\t";
        	}
        	std::cout <<std::endl;
    	}std::cout << std::endl;
}

namespace ublas = boost::numeric::ublas;

//singule value decomposition.

void householder(ublas::matrix<double> &A, ublas::matrix<double> &Q, ublas::matrix<double> &R){
    int n = A.size1();
    int m = A.size2();

    Q = ublas::identity_matrix<double>(n);
    R = A;
    for (int k = 0; k < m; ++k) {
        ublas::vector<double> x(n);
        for (int i = k; i < n; ++i) {
            x(i) = R(i,k);
        }
        double s = norm_2(x);
        if (x(k) > 0) {
            s = -s;
        }
        x(k) = x(k) - s;
        double norm_x = norm_2(x);
        if (norm_x == 0) {
            continue;
        }
        x = x / norm_x;
        ublas::matrix<double> H = ublas::identity_matrix<double>(n) - 2 * outer_prod(x,x);
        R = prod(H,R);
        Q = prod(Q,H);
    }
}

void svd_qr_shift(ublas::matrix<double> &A, ublas::matrix<double> &U, ublas::matrix<double> &S, ublas::matrix<double> &V){
    int n = A.size1();
    int m = A.size2();
    int k = std::min(n,m);
    ublas::matrix<double> Q,R;
    householder(A,Q,R);
    ublas::matrix<double> B = prod(R,Q);
    ublas::matrix<double> C = B;
    for (int i = 0; i < k; ++i) {
        double mu = C(i,i);
        ublas::matrix<double> D = C - mu * ublas::identity_matrix<double>(n);
        ublas::matrix<double> Q,R;
        householder(D,Q,R);
        C = prod(R,Q) + mu * ublas::identity_matrix<double>(n);
    }
    S = C;
    U = Q;
    V = Q;
}

int main(){
    matriz<double> a(3,3);
    a[0][0] = 0; a[0][1] = 2; a[0][2] = 3;
    a[1][0] = 4; a[1][1] = 0; a[1][2] = 0;
    a[2][0] = 7; a[2][1] = 8; a[2][2] = 9;
    printM(a);
    ublas::matrix<double> A(3,3);
    A(0,0) = 0; A(0,1) = 2; A(0,2) = 3;
    A(1,0) = 4; A(1,1) = 0; A(1,2) = 0;
    A(2,0) = 7; A(2,1) = 8; A(2,2) = 9;
    ublas::matrix<double> U(3,3), S(3,3), V(3,3);
    svd_qr_shift(A,U,S,V);
    //std::cout << U << std::endl;
    std::cout << S << std::endl;
    //std::cout << V << std::endl;

    return 0;
}