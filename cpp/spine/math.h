/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#pragma once

#include <math.h>
#ifdef MKL
	#include <mkl.h>
#endif

namespace spine
{
	namespace math
	{
		// Type definitions
		#ifdef DOUBLE_PRECISION
			typedef double realnum;
		#else
			typedef float realnum;
		#endif
		
		#ifdef MKL
			#ifdef DOUBLE_PRECISION
				typedef MKL_Complex16 complex;
			#else
				typedef MKL_Complex8 complex;
		#endif
		#else
			typedef struct {
				realnum real;
				realnum imag;
			} complex;
		#endif

		// Function prototypes - MKL optimized implementations available
		complex * allocComplex(unsigned int N);
		realnum * allocReal(unsigned int N);
		void freeComplex(complex * m);
		void freeReal(realnum * m);
		void initCopy(unsigned int N, complex * A, complex * ref);
		void add(unsigned int N, complex * a, complex * b, complex * y);
		void multiply(unsigned int dim, complex * A, complex * B, complex * C);
		void multiplyVector(unsigned int dim, complex * A, complex * b, complex * c);

		// Function prototypes
		void initZero(unsigned int N, complex * A);
		void initIdentity(unsigned int dim, complex * A);
		complex multiply(complex a, complex b);
		complex add(complex a, complex b);
		realnum abs2(complex a);
	}
}
