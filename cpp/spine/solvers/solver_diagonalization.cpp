/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "solver_diagonalization.h"

using namespace spine;

#ifdef MKL

void spine::solvers::solver_diagonalization(unsigned int dim, realnum * matrix, complex * result)
{
	// NOTE: the matrix is assumed to be real and symmetric!
	// NOTE: the matrix gets overwritten with the eigenvectors!
	unsigned int N = dim * dim;

	// 1) diagonalize the matrix: matrix = V * D * inv(V)
	realnum * D = allocReal(dim);
#ifdef DOUBLE_PRECISION
	int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', dim, matrix, dim, D);
#else
	int info = LAPACKE_ssyev(LAPACK_ROW_MAJOR, 'V', 'U', dim, matrix, dim, D);
#endif

	// 2) exponentiate, exp(-i * D)
	complex * eig = allocComplex(N);
	initZero(N, eig);
	for (unsigned int n = 0, i = 0; n < dim; n++, i += dim + 1)
	{
#ifdef DOUBLE_PRECISION
		eig[i].real = cos(D[n]);
		eig[i].imag = -sin(D[n]);
#else
		eig[i].real = cosf(D[n]);
		eig[i].imag = -sinf(D[n]);
#endif
	}

	// 3) multiply V * D * inv(V)
	complex zero = { 0.0, 0.0 };
	complex one = { 1.0, 0.0 };
	complex * tmp = allocComplex(N);
	complex * V = allocComplex(N);
	for (unsigned int n = 0; n < N; n++) {
		V[n].real = matrix[n];
		V[n].imag = 0.0;
	}

#ifdef DOUBLE_PRECISION
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &one, V, dim, eig, dim, &zero, tmp, dim);
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans, dim, dim, dim, &one, tmp, dim, V, dim, &zero, result, dim);
#else
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &one, V, dim, eig, dim, &zero, tmp, dim);
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, dim, dim, dim, &one, tmp, dim, V, dim, &zero, result, dim);
#endif

	freeReal(D);
	freeComplex(eig);
	freeComplex(tmp);
	freeComplex(V);
}


void spine::solvers::solver_diagonalization(unsigned int dim, complex * matrix, complex * result)
{
	// NOTE: the matrix is assumed to be complex and hermitian!
	// NOTE: the matrix gets overwritten with the eigenvectors!
	unsigned int N = dim * dim;

	// 1) diagonalize the matrix: matrix = V * D * inv(V)
	realnum * D = allocReal(dim);
#ifdef DOUBLE_PRECISION
	int info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', dim, matrix, dim, D);
#else
	int info = LAPACKE_cheev(LAPACK_ROW_MAJOR, 'V', 'U', dim, matrix, dim, D);
#endif

	// 2) exponentiate, exp(-i * D)
	complex * eig = allocComplex(N);
	initZero(N, eig);
	for (unsigned int n = 0, i = 0; n < dim; n++, i += dim + 1)
	{
#ifdef DOUBLE_PRECISION
		eig[i].real = cos(D[n]);
		eig[i].imag = -sin(D[n]);
#else
		eig[i].real = cosf(D[n]);
		eig[i].imag = -sinf(D[n]);
#endif
	}

	// 3) multiply V * D * inv(V)
	complex zero = { 0.0, 0.0 };
	complex one = { 1.0, 0.0 };
	complex * tmp = allocComplex(N);
	complex * V = allocComplex(N);
	for (unsigned int n = 0; n < N; n++) {
		V[n].real = matrix[n].real;
		V[n].imag = matrix[n].imag;
	}

#ifdef DOUBLE_PRECISION
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &one, V, dim, eig, dim, &zero, tmp, dim);
	for (unsigned int n = 0; n < N; n++)
		V[n].imag = -matrix[n].imag;
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans, dim, dim, dim, &one, tmp, dim, V, dim, &zero, result, dim);
#else
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &one, V, dim, eig, dim, &zero, tmp, dim);
	for (unsigned int n = 0; n < N; n++)
		V[n].imag = -matrix[n].imag;
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, dim, dim, dim, &one, tmp, dim, V, dim, &zero, result, dim);
#endif

	freeReal(D);
	freeComplex(eig);
	freeComplex(tmp);
	freeComplex(V);
}

#endif
