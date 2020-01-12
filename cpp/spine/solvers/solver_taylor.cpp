/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "solver_taylor.h"

using namespace spine;

void spine::solvers::solver_taylor(unsigned int dim, realnum * matrix, complex * result)
{
	// Default accuracy settings that can be set
	int accuracy = 5;
	int scaling = 7;

	unsigned int N = dim * dim;
	complex one = { 1.0, 0.0 };
	complex zero = { 0.0, 0.0 };

	complex * M_small = allocComplex(N);
	complex * m_exp1 = allocComplex(N);
	complex * m_exp2 = allocComplex(N);
	complex * M_power = allocComplex(N);
	complex * M_power1 = allocComplex(N);
	complex * tmpM1 = allocComplex(N);

	// Scale the matrix for more accurate calculations
	for (unsigned int i = 0; i < N; i++) {
		M_small[i].real = 0.0;

#ifdef DOUBLE_PRECISION
		M_small[i].imag = -matrix[i] / pow(2.0, (realnum)scaling);
#else
		M_small[i].imag = -matrix[i] / powf(2.0f, (realnum)scaling);
#endif

	}
	initCopy(N, M_power, M_small);

	// Exponentiation using series expansion
	initIdentity(dim, result);
	realnum factorial_i = 1.0;
	for (int i = 1; i < accuracy; i++) {
		factorial_i = factorial_i * i;

		// result = result + M_power / factorial(i);
		for (unsigned int x = 0; x < N; x++) {
			tmpM1[x].real = M_power[x].real / factorial_i;
			tmpM1[x].imag = M_power[x].imag / factorial_i;
		}
		add(N, result, tmpM1, result);

		// M_power = M_power * M_small;
		initCopy(N, M_power1, M_power);
		multiply(dim, M_power1, M_small, M_power);
	}

	// Scale back, using squaring
	for (int i = 0; i < scaling; i++) {
		initCopy(N, m_exp1, result);
		initCopy(N, m_exp2, result);
		multiply(dim, m_exp1, m_exp2, result);
	}

	freeComplex(tmpM1);
	freeComplex(M_power1);
	freeComplex(M_power);
	freeComplex(m_exp2);
	freeComplex(m_exp1);
	freeComplex(M_small);
}

void spine::solvers::solver_taylor(unsigned int dim, complex * matrix, complex * result)
{
	// Default accuracy settings that can be set
	int accuracy = 5;
	int scaling = 7;

	unsigned int N = dim * dim;
	complex one = { 1.0, 0.0 };
	complex zero = { 0.0, 0.0 };

	complex * M_small = allocComplex(N);
	complex * m_exp1 = allocComplex(N);
	complex * m_exp2 = allocComplex(N);
	complex * M_power = allocComplex(N);
	complex * M_power1 = allocComplex(N);
	complex * tmpM1 = allocComplex(N);

	// Scale the matrix for more accurate calculations
	for (unsigned int i = 0; i < N; i++) {
		
#ifdef DOUBLE_PRECISION
		M_small[i].real = matrix[i].imag / pow(2.0, (realnum)scaling);
		M_small[i].imag = -matrix[i].real / pow(2.0, (realnum)scaling);
#else
		M_small[i].real = matrix[i].imag / pow(2.0f, (realnum)scaling);
		M_small[i].imag = -matrix[i].real / powf(2.0f, (realnum)scaling);
#endif

	}
	initCopy(N, M_power, M_small);

	// Exponentiation using series expansion
	initIdentity(dim, result);
	realnum factorial_i = 1.0;
	for (int i = 1; i < accuracy; i++) {
		factorial_i = factorial_i * i;

		// result = result + M_power / factorial(i);
		for (unsigned int x = 0; x < N; x++) {
			tmpM1[x].real = M_power[x].real / factorial_i;
			tmpM1[x].imag = M_power[x].imag / factorial_i;
		}
		add(N, result, tmpM1, result);

		// M_power = M_power * M_small;
		initCopy(N, M_power1, M_power);
		multiply(dim, M_power1, M_small, M_power);
	}

	// Scale back, using squaring
	for (int i = 0; i < scaling; i++) {
		initCopy(N, m_exp1, result);
		initCopy(N, m_exp2, result);
		multiply(dim, m_exp1, m_exp2, result);
	}

	freeComplex(tmpM1);
	freeComplex(M_power1);
	freeComplex(M_power);
	freeComplex(m_exp2);
	freeComplex(m_exp1);
	freeComplex(M_small);
}
