//singule value decomposition
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

#include "mat.h"

template<class T>
inline T SQR(const T a) {return a*a;}

template<class T>
inline const T &MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

inline double MAX(const double &a, const double &b)
        {return b > a ? (b) : double(a);}

/*inline double MAX(const double &a, const double &b)
        {return b > a ? double(b) : (a);}*/

template<class T>
inline const T &MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}

inline double MIN(const double &a, const double &b)
        {return b < a ? (b) : double(a);}
/*
inline double MIN(const double &a, const double &b)
        {return b < a ? double(b) : (a);}*/

template<class T>
inline T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline double SIGN(const double &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
/*
inline double SIGN(const double &a, const double &b)
	{return (double)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}*/

template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}

//template <typename T>
class SVD{
    int m,n;
	matriz<double> u,v;
	std::vector<double> w,s;
	double eps, tsh;

    public:
	SVD(matriz<double> &a) : m(a.linhas()), n(a.colunas()), u(a), v(n,n), w(n), s(n) {
		
        eps = std::numeric_limits<double>::epsilon();
		decompose();
		reorder();
		tsh = 0.5*sqrt(m+n+1.)*w[0]*eps;
        
	}
    //refazer
	void solve(std::vector<double> &b, std::vector<double> &x,double thresh=-1.);
	//template <typename T> void solve(matriz<double> &b, std::vector<double> &x, double thresh);

	int rank(double thresh);
	int nullity(double thresh);
	matriz<double> range(double thresh);
	matriz<double> nullspace(double thresh);

	double inv_condition() {
		return (w[0] <= 0. || w[n-1] <= 0.) ? 0. : w[n-1]/w[0];
	}

	void decompose();
	void reorder();
	double pythag(const double a, const double b);

    matriz<double> getU() {return u;}
    matriz<double> getV() {return v;}
    std::vector<double> getW() {return w;}
};

void SVD::solve(std::vector<double> &b,std::vector<double> &x,double thresh) {
	int i,j,jj;
	double s;

	if (b.size() != m) throw("SVD::solve bad sizes");
	std::vector<double> tmp(n);

	tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
	for (j=0;j<n;j++) {
		s=0.0;
		if (w[j] > tsh) {
			for (i=0;i<m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=0;j<n;j++) {
		s=0.0;
		for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
}

/*template <typename T>
void SVD::solve(MatDoub_I &b, MatDoub_O &x, Doub thresh = -1.)
{
	int i,j,m=b.ncols();
	if (b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols())
		throw("SVD::solve bad sizes");
	T xx(n);
	/*for (j=0;j<m;j++) {
		for (i=0;i<n;i++) xx[i] = b[i][j];
		solve(b,xx,thresh);
		for (i=0;i<n;i++) x[i][j] = xx[i];
	}
}*/

int SVD::rank(double thresh = -1.) {
	int j,nr=0;

	tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
	for (j=0;j<n;j++) if (w[j] > tsh) nr++;

	return nr;
}

int SVD::nullity(double thresh = -1.) {
	int j,nn=0;

	tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
	for (j=0;j<n;j++) if (w[j] <= tsh) nn++;

	return nn;
}

matriz<double> SVD::range(double thresh = -1.){
	int i,j,nr=0;
	matriz<double> rnge(m,rank(thresh));

	for (j=0;j<n;j++) {
		if (w[j] > tsh) {
			for (i=0;i<m;i++) rnge[i][nr] = u[i][j];
			nr++;
		}
	}
	return rnge;
}

matriz<double> SVD::nullspace(double thresh = -1.){
	int j,jj,nn=0;
	matriz<double> nullsp(n,nullity(thresh));

	for (j=0;j<n;j++) {
		if (w[j] <= tsh) {
			for (jj=0;jj<n;jj++) nullsp[jj][nn] = v[jj][j];
			nn++;
		}
	}
	return nullsp;
}

void SVD::decompose() {
	bool flag = true;
	int i=0,its=0,j=0,jj=0,k=0,l=0,nm=0;
	double anorm=0,c=0,f=0,g=0,h=0,s=0,scale=0,x=0,y=0,z=0;

	std::vector<double> rv1(n);
	g = scale = anorm = 0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += abs(u[k][i]);
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					u[k][i] /= scale;
					s += u[k][i]*u[k][i];
				}
				f=u[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				u[i][i]=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += u[k][i]*u[k][j];
					f=s/h;
					for (k=i;k<m;k++) u[k][j] += f*u[k][i];
				}
				for (k=i;k<m;k++) u[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i+1 != n) {
			for (k=l-1;k<n;k++) scale += abs(u[i][k]);
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					u[i][k] /= scale;
					s += u[i][k]*u[i][k];
				}
				f=u[i][l-1];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				u[i][l-1]=f-g;
				for (k=l-1;k<n;k++) rv1[k]=u[i][k]/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += u[j][k]*u[i][k];
					for (k=l-1;k<n;k++) u[j][k] += s*rv1[k];
				}
				for (k=l-1;k<n;k++) u[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(abs(w[i])+abs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++)
					v[j][i]=(u[i][j]/u[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += u[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) u[i][j]=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += u[k][i]*u[k][j];
				f=(s/u[i][i])*g;
				for (k=i;k<m;k++) u[k][j] += f*u[k][i];
			}
			for (j=i;j<m;j++) u[j][i] *= g;
		} else for (j=i;j<m;j++) u[j][i]=0.0;
		++u[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=true;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if (l == 0 || abs(rv1[l]) <= eps*anorm) {
					flag=false;
					break;
				}
				if (abs(w[nm]) <= eps*anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if (abs(f) <= eps*anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=u[j][nm];
						z=u[j][i];
						u[j][nm]=y*c+z*s;
						u[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 29) throw("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=u[jj][j];
					z=u[jj][i];
					u[jj][j]=y*c+z*s;
					u[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
}

void SVD::reorder() {
	int i,j,k,s,inc=1;
	double sw;
	std::vector<double> su(m), sv(n);
	do { inc *= 3; inc++; } while (inc <= n);
	do {
		inc /= 3;
		for (i=inc;i<n;i++) {
			sw = w[i];
			for (k=0;k<m;k++) su[k] = u[k][i];
			for (k=0;k<n;k++) sv[k] = v[k][i];
			j = i;
			while (w[j-inc] < sw) {
				w[j] = w[j-inc];
				for (k=0;k<m;k++) u[k][j] = u[k][j-inc];
				for (k=0;k<n;k++) v[k][j] = v[k][j-inc];
				j -= inc;
				if (j < inc) break;
			}
			w[j] = sw;
			for (k=0;k<m;k++) u[k][j] = su[k];
			for (k=0;k<n;k++) v[k][j] = sv[k];

		}
	} while (inc > 1);
	for (k=0;k<n;k++) {
		s=0;
		for (i=0;i<m;i++) if (u[i][k] < 0.) s++;
		for (j=0;j<n;j++) if (v[j][k] < 0.) s++;
		if (s > (m+n)/2) {
			for (i=0;i<m;i++) u[i][k] = -u[i][k];
			for (j=0;j<n;j++) v[j][k] = -v[j][k];
		}
	}
}

double SVD::pythag(const double a, const double b) {
	double absa=abs(a), absb=abs(b);
	return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)) :
		(absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
}

void svd(matriz<double> &A, matriz<double> &u, matriz<double> &v, std::vector<double> &w) {
    SVD svdSolver(A);
    w = std::move(svdSolver.getW());
    u = svdSolver.getU();
    v = std::move(svdSolver.getV());

}

int main(){
    matriz <double> A(3,3,{
        1,2,3,
        4,5,6,
        7,8,9
    }),U,V;

    std::vector<double> w,s(3);
    
    svd(A,U,V,w);
    
    std::cout << std::setprecision(3) << std::fixed;
    std::cout << "A = " << A << std::endl;
    std::cout << "U = " << U << std::endl;
    for(int i = 0; i < w.size(); i++)
        std::cout << w[i] << " "; 
    std::cout << std::endl;
    std::cout << "V = " << V << std::endl;

    return 0;
}