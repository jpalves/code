#include "common.hpp"

namespace ublas = boost::numeric::ublas;

//#define DEBUG
#define CHECK_RESULT

void
householder(ublas::matrix < float >&A,
	    ublas::matrix < float >&QQ,
	    unsigned int row_start, unsigned int col_start, bool column)
{
	unsigned int size = column ? A.size1() : A.size2();
	unsigned int start = column ? row_start : col_start;

	if (start >= size)
		return;	

	ublas::vector < float > x(size);
	for (unsigned int i = 0; i < size; i++) {
		if (i < start) {
			x(i) = 0;
		} else {
			if (column)
				x(i) = A(i, col_start);
			else
				x(i) = A(row_start, i);
		}
	}

	float x_norm = norm(x);
	float alpha = sign(x(start)) * x_norm;

	// std::cout << column << " " <<  start << " x = " << x(start) << "\n";
#ifdef DEBUG
	std::cout << "||x|| = " << x_norm << "\n";
	std::cout << "alpha = " << alpha << "\n";
#endif

	ublas::vector<float> v = x;

	v(start) += alpha;
	normalize(v);

#ifdef DEBUG
	std::cout << "v = " << v << "\n";
#endif

	ublas::matrix < float > Q;

#ifdef DEBUG
	Q.resize(size, size);
	eye(Q);

	for (unsigned int i = start; i < Q.size1(); i++) {
		for (unsigned int j = start; j < Q.size2(); j++) {
			Q(i, j) = Q(i, j) - 2 * v(i) * v(j);
		}
	}

	std::cout << "Q  = " << Q << "\n";
#endif

	if (column) {
		for (unsigned int i = 0; i < A.size2(); i++) {
			float sum_Av = 0.0f;

			for (unsigned int j = 0; j < A.size1(); j++)
				sum_Av = sum_Av + (v(j) * A(j, i));
			for (unsigned int j = 0; j < A.size1(); j++)
				A(j, i) = A(j, i) - 2 * v(j) * sum_Av;
		}

		for (unsigned int i = 0; i < A.size1(); i++) {
			float sum_Qv = 0.0f;

			for (unsigned int j = 0; j < A.size1(); j++)
				sum_Qv = sum_Qv + (v(j) * QQ(i, j));
			for (unsigned int j = 0; j < A.size1(); j++)
				QQ(i, j) = QQ(i, j) - 2 * v(j) * sum_Qv;
		}

	} else {
		for (unsigned int i = 0; i < A.size1(); i++) {
			float sum_Av = 0.0f;

			for (unsigned int j = 0; j < A.size2(); j++)
				sum_Av = sum_Av + (v(j) * A(i, j));
			for (unsigned int j = 0; j < A.size2(); j++)
				A(i, j) = A(i, j) - 2 * v(j) * sum_Av;
		}

		for (unsigned int i = 0; i < A.size2(); i++) {
			float sum_Qv = 0.0f;

			for (unsigned int j = 0; j < A.size2(); j++)
				sum_Qv = sum_Qv + (v(j) * QQ(i, j));
			for (unsigned int j = 0; j < A.size2(); j++)
				QQ(i, j) = QQ(i, j) - 2 * v(j) * sum_Qv;
		}
	}
}

void
svd_qr_shift(ublas::matrix < float >&u,
	     ublas::matrix < float >&v,
	     ublas::vector < float >&q, ublas::vector < float >&e)
{
	int n = q.size();
	int m = u.size1();

	std::cout << u.size1() << " " << u.size2() << "\n";

	bool goto_test_conv = false;

	for (int k = n - 1; k >= 0; k--) {
		//std::cout << "U = " << u << std::endl;

		for (int iter = 0; iter < ITER_MAX; iter++) {
			// test for split
			int l;
			for (l = k; k >= 0; l--) {
				goto_test_conv = false;
				if (fabs(e[l]) <= EPS) {
					// set it
					goto_test_conv = true;
					break;
				}

				if (fabs(q[l - 1]) <= EPS) {
					// goto
					break;
				}
			}

			if (!goto_test_conv) {
				float c = 0.0;
				float s = 1.0;

				int l1 = l - 1;

				for (int i = l; i <= k; i++) {
					float f = s * e[i];
					e[i] = c * e[i];

					if (fabs(f) <= EPS) {
						break;
					}

					float g = q[i];
					float h = pythag(f, g);
					q[i] = h;
					c = g / h;
					s = -f / h;

					for (int j = 0; j < m; j++) {
						float y = u(j, l1);
						float z = u(j, i);
						u(j, l1) = y * c + z * s;
						u(j, i) = -y * s + z * c;
					}
				}
			}

			float z = q[k];

			if (l == k) {
				if (z < 0.0f) {
					q[k] = -z;

					for(int j = 0; j < n; j++)
						v(j, k) = -v(j,k);
				}

				break;
			}

			if (iter >= ITER_MAX - 1) {
				break;
			}

			float x = q[l];
			float y = q[k - 1];
			float g = e[k - 1];
			float h = e[k];
			float f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			
			g = pythag(f, 1.0);

			if (f < 0) {
				f = ((x - z) * (x + z) + h * (y / (f - g) - h)) / x;
			} else {
				f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
			}

			float c = 1.0;
			float s = 1.0;

			for (int i = l + 1; i <= k; i++) {
				g = e[i];
				y = q[i];
				h = s * g;
				g = c * g;
				float z = pythag(f, h);
				e[i - 1] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = -x * s + g * c;
				h = y * s;
				y = y * c;

				for (int j = 0; j < n; j++) {
					x = v(j, i - 1);
					z = v(j, i);
					v(j, i - 1) = x * c + z * s;
					v(j, i) = -x * s + z * c;
				}

				z = pythag(f, h);
				q[i - 1] = z;
				c = f / z;
				s = h / z;
				f = c * g + s * y;
				x = -s * g + c * y;

				for (unsigned int j = 0; j < m; j++) {
					float y = u(j, i - 1);
					float z = u(j, i);
					u(j, i - 1) = y * c + z * s;
					u(j, i) = -y * s + z * c;
				}
			}
			e[l] = 0.0;
			e[k] = f;
			q[k] = x;
		}
	}
}

void
svd(ublas::matrix < float >&A,
    	ublas::matrix < float >&QQL,
    	ublas::matrix < float >&QQW, 
		ublas::matrix < float >&QQR)
{
	unsigned int row_num = A.size1();
	unsigned int col_num = A.size2();

	QQL.resize(row_num, row_num);
    QQW.resize(row_num, col_num);
	QQR.resize(col_num, col_num);

	eye(QQL);
	eye(QQR);

	unsigned int to = std::min(row_num, col_num);

	for (unsigned int i = 0; i < to; i++) {
		householder(A, QQL, i, i, true);
		householder(A, QQR, i, i + 1, false);
	}

	ublas::vector < float > d(to, 0.0f);
	ublas::vector < float > s(to + 1, 0.0f);

	for (unsigned int i = 0; i < to; i++) {
		d(i) = A(i, i);
		if (i < (to - 1))
			s(i + 1) = A(i, i + 1);
	}

	// std::cout << d << "\n";
	// std::cout << s << "\n";

	// std::cout << QQL << std::endl;
	// std::cout << QQR << std::endl;

	svd_qr_shift(QQL, QQR, d, s);

	std::cout << "Sigmas: " << d << "\n";

	QQW.clear();

	for (unsigned int i = 0; i < to; i++)
		QQW(i, i) = d(i);

}

bool check_bidiag(ublas::matrix < float >&A)
{
	const float EPS = 0.0001f;

	for (unsigned int i = 0; i < A.size1(); i++) {
		for (unsigned int j = 0; j < A.size2(); j++) {
			if ((std::abs(A(i, j)) > EPS) && (i != j)
			    && ((i + 1) != j)) {
				std::
				    cout << "Failed at " << i << " " << j << " "
				    << A(i, j)
				    << "\n";
				return false;
			}
		}
	}
	return true;
}

int main()
{
	ublas::matrix < float > in;

	std::fstream f;
	f.open("data/wiki.example", std::fstream::in);
	f >> in;
	f.close();

//  srand((unsigned int)time(0));
//	random_fill(in, 250, 150);

	ublas::matrix < float > ref = in;

	pretty_print("Input:", in);

#ifdef DEBUG
	std::cout << in << "\n";
#endif
	ublas::matrix < float > QQL;
    ublas::matrix < float > QQW;
	ublas::matrix < float > QQR;

    svd(in, QQL, QQW, QQR);

	// std::cout << in << "\n";
	ublas::matrix < float > result;

#ifdef CHECK_RESULT
	// std::cout << "QL: " << QQL << "\n";
 //    std::cout << "QW: " << QQW << "\n";
	// std::cout << "QR: " << QQR << "\n";
	// std::cout << "A:: " << in << "\n";

    result = ublas::prod(QQW, trans(QQR));
	result = ublas::prod(QQL, result);

	// std::cout << result << "\n";
#endif

	pretty_print("Bidiag:", in);


	std::cout << "DIFF    = " << matrix_compare(result, ref) << "\n";
	std::cout << "Is bidiag " << check_bidiag(in) << "\n";

	return 0;
}
