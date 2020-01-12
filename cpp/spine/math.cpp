/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "math.h"

using namespace spine::math;

#ifdef MKL

complex * spine::math::allocComplex(unsigned int N)
{
	return (complex *)mkl_malloc(N * sizeof(complex), 64);
}

realnum * spine::math::allocReal(unsigned int N)
{
	return (realnum *)mkl_malloc(N * sizeof(realnum), 64);
}

void spine::math::freeComplex(complex * m)
{
	mkl_free(m);
}

void spine::math::freeReal(realnum * m)
{
	mkl_free(m);
}

void spine::math::initCopy(unsigned int N, complex * A, complex * ref)
{
#ifdef DOUBLE_PRECISION
	cblas_zcopy(N, ref, 1, A, 1);
#else
	cblas_ccopy(N, ref, 1, A, 1);
#endif
}

void spine::math::add(unsigned int N, complex * a, complex * b, complex * y)
{
#ifdef DOUBLE_PRECISION
	vzAdd(N, a, b, y);
#else
	vcAdd(N, a, b, y);
#endif
}

void spine::math::multiply(unsigned int dim, complex * A, complex * B, complex * C)
{
	complex zero = { 0.0, 0.0 };
	complex one = { 1.0, 0.0 };

	// C = A * B
#ifdef DOUBLE_PRECISION
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &one, A, dim, B, dim, &zero, C, dim);
#else
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &one, A, dim, B, dim, &zero, C, dim);
#endif
}

void spine::math::multiplyVector(unsigned int dim, complex * A, complex * b, complex * c)
{
	complex zero = { 0.0, 0.0 };
	complex one = { 1.0, 0.0 };

	// c = A * b
#ifdef DOUBLE_PRECISION
	cblas_zgemv(CblasRowMajor, CblasNoTrans, dim, dim, &one, A, dim, b, 1, &zero, c, 1);
#else
	cblas_cgemv(CblasRowMajor, CblasNoTrans, dim, dim, &one, A, dim, b, 1, &zero, c, 1);
#endif
}

#else

complex * spine::math::allocComplex(unsigned int N)
{
	return new complex[N];
}

realnum * spine::math::allocReal(unsigned int N)
{
	return new realnum[N];
}

void spine::math::freeComplex(complex * m)
{
	delete[] m;
}

void spine::math::freeReal(realnum * m)
{
	delete[] m;
}

void spine::math::initCopy(unsigned int N, complex * A, complex * ref)
{
	for (unsigned int n = 0; n < N; n++) {
		A[n].real = ref[n].real;
		A[n].imag = ref[n].imag;
	}
}

void spine::math::add(unsigned int N, complex * a, complex * b, complex * y)
{
	for (unsigned int n = 0; n < N; n++) {
		y[n].real = a[n].real + b[n].real;
		y[n].imag = a[n].imag + b[n].imag;
	}
}

void spine::math::multiply(unsigned int dim, complex * A, complex * B, complex * C)
{
	// C = A * B
	for (unsigned int i = 0; i < dim; i++) {
		for (unsigned int j = 0; j < dim; j++) {
			C[j + dim*i].real = 0;
			C[j + dim*i].imag = 0;
			for (unsigned int k = 0; k < dim; k++) {
				C[j + dim*i].real += A[k + dim*i].real * B[j + dim*k].real - A[k + dim*i].imag * B[j + dim*k].imag;
				C[j + dim*i].imag += A[k + dim*i].real * B[j + dim*k].imag + A[k + dim*i].imag * B[j + dim*k].real;
			}
		}
	}
}

void spine::math::multiplyVector(unsigned int dim, complex * A, complex * b, complex * c)
{
	// c = A * b
	for (unsigned int i = 0; i < dim; i++) {
		c[i].real = 0;
		c[i].imag = 0;
		for (unsigned int k = 0; k < dim; k++) {
			c[i].real += A[k + dim * i].real * b[k].real - A[k + dim * i].imag * b[k].imag;
			c[i].imag += A[k + dim * i].real * b[k].imag + A[k + dim * i].imag * b[k].real;
		}
	}
}

#endif

void spine::math::initZero(unsigned int N, complex * A)
{
	for (unsigned int n = 0; n < N; n++) {
		A[n].real = (realnum) 0.0;
		A[n].imag = (realnum) 0.0;
	}
}

void spine::math::initIdentity(unsigned int dim, complex * A)
{
	unsigned int N = dim * dim;

	initZero(N, A);
	for (unsigned int n = 0; n < N; n += dim + 1)
		A[n].real = (realnum) 1.0;
}

complex spine::math::multiply(complex a, complex b)
{
	complex result;
	result.real = a.real * b.real - a.imag * b.imag;
	result.imag = a.real * b.imag + a.imag * b.real;
	return result;
}

complex spine::math::add(complex a, complex b)
{
	complex result;
	result.real = a.real + b.real;
	result.imag = a.imag + b.imag;
	return result;
}

realnum spine::math::abs2(complex a)
{
	return a.real * a.real + a.imag * a.imag;
}
