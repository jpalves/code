#include <iostream>
#include <cmath>
#include "mat.h"
#include <iomanip>
#include <sstream>

#define sign(x)  0.0 ? 1.0 : -1.0
static const int ITER_MAX = 50;
static const double EPS = 0.00001;

float norm(matriz<double>&x)
{
    float x_norm = 0.0;
    for (unsigned int i = 0; i < x.linhas(); i++)
        x_norm += std::pow(x[i][0], 2);
    x_norm = std::sqrt(x_norm);
    return x_norm;
}

void pretty_print(const char* tag, matriz<double>& m)
{
    printf("%s\n", tag);
    for(std::size_t i = 0; i < m.linhas(); i++)
    {
        printf("\t");
        for(std::size_t j = 0; j < m.colunas(); j++)
        {
            printf("%5.2f\t", (fabs(m[i][j]) < 0.00001f?0.0f:m[i][j]));
        }
        printf("\n");
    }
}
/*
size_t number_of_digits(double n=5) {
	std::ostringstream strs;

	strs << n;
	return strs.str().size();
}*/

void print_matrix(matriz <double> M) {
	size_t n = M.colunas(), m = M.linhas();
	int nmax= 6;
	size_t max_len_per_column[nmax];

	for (size_t j = 0; j < m; ++j) {
		size_t max_len {};

		for (size_t i = 0; i < n; ++i)
			if (const auto num_length {number_of_digits(M[i][j])}; num_length > max_len)
				max_len = num_length;

		max_len_per_column[j] = max_len;
	}

	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			std::cout << (j == 0 ? "\n| " : "") << std::setw(max_len_per_column[j]) << (fabs(M[i][j]) < 0.00001f?0.0f:M[i][j]) << (j == m - 1 ? " |" : " ");

	std::cout << '\n';
}

void normalize(matriz<double> &x){
    float x_norm = norm(x);

    for (unsigned int i = 0; i < x.linhas(); i++) {
        x[i][0] /= x_norm;
    }
}

double pythag(double a, double b){
    double absa = fabs(a);
    double absb = fabs(b);

    if (absa > absb) {
        return absa * sqrt(1.0f + pow(absb / absa, 2));
    } else {
        return absb * sqrt(1.0f + pow(absa / absb, 2));
    }
}

//#define DEBUG
#define CHECK_RESULT

void householder(matriz <double> &A,
	             matriz <double> &QQ,
	             unsigned int row_start, unsigned int col_start, bool column){
	unsigned int size = column ? A.linhas() : A.colunas();
	unsigned int start = column ? row_start : col_start;

	if (start >= size)
		return;	

	matriz <double> v(size,1,0.0);
	for (unsigned int i = 0; i < size; i++) {
		if (i < start) {
			v.asVector()[i] = 0;
		} else {
			if (column)
				v.asVector()[i] = A[i][col_start];
			else
				v.asVector()[i] = A[row_start][i];
		}
	}
    //std::cout << x << "\n";

	float x_norm = norm(v);
	float alpha  = sign(v.asVector()[start]) * x_norm;

	// std::cout << column << " " <<  start << " x = " << x(start) << "\n";

	//matriz<double> v = x;

	//v.set 

	v.asVector()[start] += alpha;
	normalize(v);

	//matriz <double> Q;


	if (column) {
		for (unsigned int i = 0; i < A.colunas(); i++) {
			double sum_Av = 0.0f;

			for (unsigned int j = 0; j < A.linhas(); j++)
				sum_Av = sum_Av + (v.asVector()[j] * A[j][i]);
			for (unsigned int j = 0; j < A.linhas(); j++)
				A[j][i] = A[j][i] - 2 * v.asVector()[j] * sum_Av;
		}

		for (unsigned int i = 0; i < A.linhas(); i++) {
			double sum_Qv = 0.0f;

			for (unsigned int j = 0; j < A.linhas(); j++)
				sum_Qv = sum_Qv + (v.asVector()[j] * QQ[i][j]);
			for (unsigned int j = 0; j < A.linhas(); j++)
				QQ[i][j] = QQ[i][j] - 2 * v.asVector()[j] * sum_Qv;
		}

	} else { //estou aqui
		//pretty_print("debug:", A);
		for (unsigned int i = 0; i < A.linhas(); i++) {
			double sum_Av = 0.0f;

			for (unsigned int j = 0; j < A.colunas(); j++)
				sum_Av = sum_Av + (v[j][0] * A[i][j]);
			for (unsigned int j = 0; j < A.colunas(); j++)
				A[i][j] = A[i][j] - 2 * v[j][0] * sum_Av;
		}

		for (unsigned int i = 0; i < A.colunas(); i++) {
			double sum_Qv = 0.0f;

			for (unsigned int j = 0; j < A.colunas(); j++)
				sum_Qv = sum_Qv + (v[j][0] * QQ[i][j]);
			for (unsigned int j = 0; j < A.colunas(); j++)
				QQ[i][j] = QQ[i][j] - 2 * v[j][0] * sum_Qv;
		}
	}
}

void svd_qr_shift(matriz<double>&u,
	     matriz<double>&v,
	     matriz<double>&q, matriz<double>&e){
	int n = q.linhas();
	int m = u.linhas();

	//std::cout << u.linhas() << " " << u.colunas() << "\n";

	bool goto_test_conv = false;

	for (int k = n - 1; k >= 0; k--) {
		//std::cout << "U = " << u << std::endl;

		for (int iter = 0; iter < ITER_MAX; iter++) {
			// test for split
			int l;
			for (l = k; k >= 0; l--) {
				goto_test_conv = false;
				if (fabs(e[l][0]) <= EPS) {
					// set it
					goto_test_conv = true;
					break;
				}

				if (fabs(q[l - 1][0]) <= EPS) {
					// goto
					break;
				}
			}

			if (!goto_test_conv) {
				double c = 0.0;
				double s = 1.0;

				int l1 = l - 1;

				for (int i = l; i <= k; i++) {
					double f = s * e[i][0];
					e[i][0] = c * e[i][0];

					if (fabs(f) <= EPS) {
						break;
					}

					double g = q[i][0];
					double h = pythag(f, g);
					q[i][0] = h;
					c = g / h;
					s = -f / h;

					for (int j = 0; j < m; j++) {
						double y = u[j][l1];
						double z = u[j][i];
						u[j][l1] = y * c + z * s;
						u[j][i] = -y * s + z * c;
					}
				}
			}

			double z = q[k][0];

			if (l == k) {
				if (z < 0.0f) {
					q[k][0] = -z;

					for(int j = 0; j < n; j++)
						v[j][k] = -v[j][k];
				}

				break;
			}

			if (iter >= ITER_MAX - 1) {
				break;
			}

			double x = q[l][0];
			double y = q[k - 1][0];
			double g = e[k - 1][0];
			double h = e[k][0];
			double f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			
			g = pythag(f, 1.0);

			if (f < 0) {
				f = ((x - z) * (x + z) + h * (y / (f - g) - h)) / x;
			} else {
				f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
			}

			double c = 1.0;
			double s = 1.0;

			for (int i = l + 1; i <= k; i++) {
				g = e[i][0];
				y = q[i][0];
				h = s * g;
				g = c * g;
				double z = pythag(f, h);
				e[i - 1][0] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = -x * s + g * c;
				h = y * s;
				y = y * c;

				for (int j = 0; j < n; j++) {
					x = v[j][i - 1];
					z = v[j][i];
					v[j][i - 1] = x * c + z * s;
					v[j][i] = -x * s + z * c;
				}

				z = pythag(f, h);
				q[i - 1][0] = z;
				c = f / z;
				s = h / z;
				f = c * g + s * y;
				x = -s * g + c * y;

				for (unsigned int j = 0; j < m; j++) {
					double y = u[j][i - 1];
					double z = u[j][i];
					u[j][i - 1] = y * c + z * s;
					u[j][i] = -y * s + z * c;
				}
			}
			e[l][0] = 0.0;
			e[k][0] = f;
			q[k][0] = x;
		}
	}
}
// ---------------------------------------------------------------------------
void svd(matriz<double>&A,
    	matriz<double>&QQL,
    	matriz<double>&QQW, 
		matriz<double>&QQR){
	int row_num = A.linhas();
	int col_num = A.colunas();

	QQL = zeros<double>(row_num, row_num);
    //QQW = zeros<double>(row_num, col_num);
    QQR = zeros<double>(col_num, col_num);

    //QQW.resize(row_num, col_num);
	//QQR.resize(col_num, col_num);

	QQL = QQL.eye();
	QQR = QQR.eye();
    std::cout << "linhas: " << row_num << std::endl;
	std::cout << "coluna: " << col_num << std::endl;
	int to = std::min(row_num, col_num);

	for (int i = 0; i < to; i++) {
		householder(A, QQL, i, i, true);
		householder(A, QQR, i, i + 1, false);
	}

	matriz <double> d(to, 1,0.0f);
	matriz <double> s(4, 1,0.0f);

	for (int i = 0; i < to; i++) {
		d[i][0] = A[i][i];
		if (i < (to - 1))
			s[i + 1][0] = A[i][i + 1];
	}

	//std::cout <<to<<std::endl<< d << "\n";
	//std::cout << s << "\n";
	//std::cout << QQL << std::endl;
	//std::cout << QQR << std::endl;

	svd_qr_shift(QQL, QQR, d, s);

	std::cout << "Sigmas: " << d << "\n";
    
	//QQW.clear(); // = zeros(row_num, col_num); modificar
    //QQW = zeros<double>(row_num, col_num);
	//for (int i = 0; i < to; i++)
	QQW = std::move(d);

}

bool check_bidiag(matriz <double>&A)
{
	const float EPS = 0.0001f;

	for (unsigned int i = 0; i < A.linhas(); i++) {
		for (unsigned int j = 0; j < A.colunas(); j++) {
			if ((std::abs(A[i][j]) > EPS) && (i != j)
			    && ((i + 1) != j)) {
				std::
				    cout << "Failed at " << i << " " << j << " "
				    << A[i][j]
				    << "\n";
				return false;
			}
		}
	}
	return true;
}



int main(){
	matriz<double> in(3, 3,{1,2,3,4,5,6,7,8,9});
	//ublas::matrix < float > ref = in;

	pretty_print("Input:", in);

	matriz <double> QQL; //U
    matriz <double> QQW; //s 
	matriz <double> QQR; //V

    svd(in, QQL, QQW, QQR);

	// std::cout << in << "\n";
	matriz <double> result;



	pretty_print("Bidiag:", in);
    pretty_print("QQL:", QQL);
    pretty_print("QQW:", QQW);
    pretty_print("QQR:", QQR);

	print_matrix(in);
	//std::cout << "DIFF    = " << matrix_compare(result, ref) << "\n";
	///std::cout << "Is bidiag " << check_bidiag(in) << "\n";

	return 0;
}
