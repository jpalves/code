#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include "mat.h"


using namespace std;

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

void QR_factorization(matriz<complexo<double>> A, matriz<complexo<double>> &Q,matriz<complexo<double>> &R){
    int m = A.linhas();
    int n = A.colunas();

    Q =  zeros<complexo<double>>(m,n);
    R =  zeros<complexo<double>>(n,n);

    for(int i = 0; i < n; i++){
        matriz<complexo<double>> v = A.coluna(i);
    
        for(int j = 0; j < i; j++){
            matriz<complexo<double>> q = Q.coluna(j);
            
            R[j][i] = (q.transposta()*v)[0][0];
            v = v - R[j][i]*q;
        }
        R[i][i] = v.norma();
        v = v/R[i][i];
        for (int k = 0; k < n; k++){
            Q[k][i] = v.asVector()[k];
        }
    }
}

int main(){
    matriz<complexo<double>> A(3,3,{0 , 2 , 3 ,
                                    4 , 0 , 0 , 
                                    7 , 8 , 9});
    matriz<complexo<double>> Q(3,3);
    matriz<complexo<double>> R(3,3);

    QR_factorization(A,Q,R);
    cout << Q << endl;
    cout << R << endl;
    
    for (int i = 0; i < 20; i++){
        QR_factorization(A,Q,R);
        A = R*Q;
    }
    cout << A << endl;

    return 0;
}