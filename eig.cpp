#include <iostream>
#include <iomanip>
#include <eigen3/Eigen/Dense>
#include <complex>
#include <mat.h>

#if (defined __GNUC__) && (__GNUC__>4 || __GNUC_MINOR__>=7)
  #undef _GLIBCXX_ATOMIC_BUILTINS
  #undef _GLIBCXX_USE_INT128
#endif

template <typename T> class complexo:public std::complex<T>{ 
private:
   template<typename U> friend std::ostream &operator<<(std::ostream &stream,complexo<U> z);
   template<typename U> friend complexo operator *(double real,complexo<U> z){
        return complexo<U>(real*z.real(),real*z.imag());
   }
   template<typename U>friend complexo operator *(complexo<U> z,double real){
        return complexo<U>(real*z.real(),real*z.imag());
   }
   template<typename U> friend complexo operator +(double real,complexo<U> z){
          return complexo<U>(real+z.real(),z.imag());
   }       
   template<typename U> friend complexo operator +(complexo<U> z,double real){
          return complexo<U>(real+z.real(),z.imag());
   }
   template<typename U> friend complexo operator -(double real,complexo<U> z){
          return complexo<U>(real-z.real(),-z.imag());
   }       
   template<typename U> friend complexo operator -(complexo<U> z,double real){
          return complexo<U>(-real+z.real(),z.imag());
   }
   /*
   template<typename U> friend complexo operator /(double real,complexo<U> z){
             return complexo<U>(real/z.real(),real/z.imag());
   }
   template<typename U> friend complexo operator /(complexo<U> z,double real){
                 return complexo<U>(z.real()/real,z.imag()/real);
   }*/
public:
    complexo(T r=0,T i=0):std::complex<T>(r,i){}
    complexo(const std::complex<T> &c):std::complex<T>(c){}
    /*complexo<T> operator*(const complexo<T> &c){
        return complexo<T>(this->real()*c.real()-this->imag()*c.imag(),this->real()*c.imag()+this->imag()*c.real());
    }*//*
    complexo<T> operator+(const complexo<T> &c){
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
    /*
    complexo<T> operator/(const complexo<T> &c){
        return complexo<T>((this->real()*c.real()+this->imag()*c.imag())/(c.real()*c.real()+c.imag()*c.imag()),(this->imag()*c.real()-this->real()*c.imag())/(c.real()*c.real()+c.imag()*c.imag()));
    }*/
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
    stream <<std::setprecision(5)<<std::setiosflags(std::ios::fixed);   		
	if(!(z.imag()-1)&&!z.real()){ stream << "j"; return stream;}
	if(!(z.imag()+1)&&!z.real()){ stream <<"-j"; return stream;}
	if(!(z.imag()-1))           { stream <<z.real()<<"+j"; return stream;}
	if(!(z.imag()+1))           { stream <<z.real()<<"-j"; return stream;}
	//if(!z.imag())               { stream << z.real();return stream;}
    //if(!z.real())               { stream <<z.imag()<<"j"; return stream;}
	if(z.imag() < 0) stream <<z.real()<<z.imag()<<'j';
	else             stream <<z.real()<<"+"<<z.imag()<<'j';
	stream <<std::resetiosflags(std::ios::fixed);

    return stream;
}



int main(){
    Eigen::Matrix3d A;
    A << 0, 2, 3,
         4, 0, 0,
         7, 8, 9;

    matriz <complexo<double>> D(3,3),V(3,3),B(3,3);
    matriz <double> LOG(3,3);
    Eigen::EigenSolver<Eigen::Matrix3d> es(A);
    
    for (int i = 0; i < D.linhas(); i++){
        D[i][i] = log(es.eigenvalues()(i));
    }
    
    for (int i = 0; i < V.linhas(); i++){
        for (int j = 0; j < V.colunas(); j++){
            V[i][j] = es.eigenvectors()(i,j);
        }
    }
    
    std::cout << D << std::endl;
    std::cout << V << std::endl;
    /*matriz <double> C = 
    s*/

    B = V * D * V.inv();

    for (int i = 0; i < V.linhas(); i++){
        for (int j = 0; j < V.colunas(); j++)
            LOG[i][j] = B[i][j].real();
    }


    std::cout << LOG << std::endl;
    return 0;
}