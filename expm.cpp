//exponential matrix.

#include <iostream>
#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>
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

template <typename T> class complexo:public std::complex<T>{ 
private:
   template<typename U> friend std::ostream &operator<<(std::ostream &stream,complexo<U> z);
   template<typename U> friend complexo operator *(double real,complexo<U> z){
        return complexo(real*z.real(),real*z.imag());
   }
   template<typename U>friend complexo operator *(complexo<U> z,double real){
        return complexo(real*z.real(),real*z.imag());
   }
   template<typename U> friend complexo operator +(double real,complexo<U> z){
          return complexo(real+z.real(),z.imag());
   }       
   template<typename U> friend complexo operator +(complexo<U> z,double real){
          return complexo(real+z.real(),z.imag());
}
public:
    complexo(T r=0,T i=0):std::complex<T>(r,i){}
    complexo(const std::complex<T> &c):std::complex<T>(c){}
    /*complexo<T> operator*(const complexo<T> &c){
        return complexo<T>(this->real()*c.real()-this->imag()*c.imag(),this->real()*c.imag()+this->imag()*c.real());
    }*/
    /*complexo<T> operator+(const complexo<T> &c){
        return complexo<T>(this->real()+c.real(),this->imag()+c.imag());
    }*/
    complexo<T> operator-(const complexo<T> &c){
        return complexo<T>(this->real()-c.real(),this->imag()-c.imag());
    }
    complexo<T> operator*(const T s){
        return complexo<T>(this->real()*s,this->imag()*s);
    }
    complexo<T> operator/(const T s){
        return complexo<T>(this->real()/s,this->imag()/s);
    }
    complexo<T> operator/(const complexo<T> &c){
        return complexo<T>((this->real()*c.real()+this->imag()*c.imag())/(c.real()*c.real()+c.imag()*c.imag()),(this->imag()*c.real()-this->real()*c.imag())/(c.real()*c.real()+c.imag()*c.imag()));
    }
    /*complexo operator  *(complexo<T> z){
        return complexo<T>(this->real()*z.real()-this->imag()*z.imag(),this->real()*z.imag()+this->imag()*z.real());
    }*/
    /*
    complexo<T> operator^(const int n){
        complexo<T> c(1,0);
        for (int i = 0; i < n; ++i) {
            c = c*(*this);
        }
        return c;
    }*/
    
};

template <typename U> std::ostream &operator <<(std::ostream &stream,complexo<U> z){
        		
	if(!(z.imag()-1)&&!z.real()){ stream << 'j'; return stream;}
	if(!(z.imag()+1)&&!z.real()){ stream <<"-j"; return stream;}
	if(!(z.imag()-1))           { stream <<z.real()<<"+j"; return stream;}
	if(!(z.imag()+1))           { stream <<z.real()<<"-j"; return stream;}
	if(!z.imag())               { stream << z.real(); return stream;}
	if(!z.real())               { stream <<z.imag()<<"j"; return stream;}
	if(z.imag() < 0) stream << z.real()<<z.imag()<<'j';
	else             stream << z.real()<<'+'<<z.imag()<<'j';
	
    return stream;
}

complexo<double> fact(complexo<double> n){
	return (n != 0.0?n*fact(n-1):1);
}
/*
template <typename T>matriz<T> zeros(int linhas,int colunas){
    matriz <T> temp(linhas,colunas);
    
    for(int i=0;i < temp.linhas();i++)
        for(int j=0;j < temp.colunas();j++)
            temp[i][j] = 0;
    return temp;
}*/

template <typename T>matriz<T> eye(int dim){
	matriz <T> temp = zeros<T>(dim,dim);
	
	for(int i=0;i < temp.linhas();i++)
		temp[i][i] = 1;
	return temp;
}

matriz<complexo<double>> powm(matriz<complexo<double>> A,int e){
	matriz<complexo<double>> temp = eye<complexo<double>>(A.linhas());
 	
	for(int i=e;i > 0;i--)
		temp = temp*A;
	return temp;
}

matriz<complexo<double>> expm(matriz<complexo<double>> A){
	matriz <complexo<double>> temp(A.linhas(),A.colunas(),complexo<double>(0));
	
	for(int i=0;i < 100;i++)
		temp = temp + powm(A,i)/fact(i);
	return temp;					
}

template <typename T> matriz <T> diag(std::initializer_list<T> in){
    matriz <T> temp(in.size(),in.size());
    int i=0;
    for(auto it = in.begin();it != in.end();it++,i++)
        temp[i][i] = *it;
    return temp;
}

//corrigir
template <typename T> matriz <T> diag(std::vector <T> v){
    matriz <T> temp(v.size(),v.size());
    for(int i=0;i < v.size();i++)
        temp[i][i] = v[i];
    return temp;
}

void spdiags(std::vector<double> &d, std::vector<int> &k, matriz <complexo<double>> &A){
    int n = A.linhas();
    int m = A.colunas();
    int p = d.size();
    int q = k.size();
    if (p != q){
        std::cout << "The size of d and k must be the same." << std::endl;
        return;
    }
    for (int i = 0; i < p; i++){
        int j = k[i];
        if (j >= 0){
            for (int l = 0; l < n - j; l++){
                A[l][l + j] = d[i];
            }
        }
        else{
            for (int l = 0; l < m + j; l++){
                A[l - j][l] = d[i];
            }
        }
    }
}

void spdiags(std::initializer_list<double> d, std::initializer_list<double> k, matriz <complexo<double>> &A){
    int n = A.linhas();
    int m = A.colunas();
    int p = d.size();
    int q = k.size();
    if (p != q){
        std::cout << "The size of d and k must be the same." << std::endl;
        return;
    }
    int i = 0;
    for (auto it = k.begin(); it != k.end(); it++, i++){
        int j = *it;
        if (j >= 0){
            for (int l = 0; l < n - j; l++){
                A[l][l + j] = d.begin()[i];
            }
        }
        else{
            for (int l = 0; l < m + j; l++){
                A[l - j][l] = d.begin()[i];
            }
        }
    }
}

using namespace std;

int main(){
    complexo <double> a,j = sqrt(complexo<double>(-1));
    matriz <complexo<double>> A(3,3,{
        0,2,3,
        4,5,6,
        7,8,9}
    ),B=A.eye(),C,D(3,3,5);

    C = diag<complexo<double>>({1,2,3});
    printM(expm(A));

    a = 5 + 2*j; 
    spdiags({ -1, 2, -1 }, { -1, 0, 1 }, D);

    cout << fact(a.real())<<" -> "<<a<<endl;
    cout << endl;
    //std::cout << std::setprecision(3) << std::fixed;
    printM(C);
    printM(B);
    printM(D);

    complexo<double> *k = B.asVector();
    k[2] = 15;

    for(int i=0;i < B.linhas()*B.colunas();i++)
        cout << k[i] << " ";
    cout << endl;
    //printM(treta);
    cout << A.inv();
    return 0;
}
