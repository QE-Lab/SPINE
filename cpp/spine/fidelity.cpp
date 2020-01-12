/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "fidelity.h"

using namespace spine;

realnum spine::fidelity(unsigned int dim, complex * U, complex * Uideal)
{
	// Formula below is only valid for 2x2 operations
	assert(dim == 2);

	// Return the process fidelity
	complex F = { 0, 0 };
	for (unsigned int i = 0; i < 4; i++) {
		Uideal[i].imag = -Uideal[i].imag;
		F = add(F, multiply(Uideal[i], U[i]));
	}
	return abs2(F) / ((realnum) 4.0);
}

realnum spine::fidelity(unsigned int dim, complex * U, realnum theta, realnum phi)
{
	// Calculate Uideal here
	complex Uideal[4];

#ifdef DOUBLE_PRECISION
	Uideal[0].real = cos(theta / 2.0);
	Uideal[0].imag = 0;
	Uideal[1].real = sin(phi) * sin(theta / 2.0);
	Uideal[1].imag = -cos(phi) * sin(theta / 2.0);
#else
	Uideal[0].real = cosf(theta / 2.0f);
	Uideal[0].imag = 0;
	Uideal[1].real = sinf(phi) * sinf(theta / 2.0f);
	Uideal[1].imag = -cosf(phi) * sinf(theta / 2.0f);
#endif
	Uideal[2].real = -Uideal[1].real;
	Uideal[2].imag = Uideal[1].imag;
	Uideal[3].real = Uideal[0].real;
	Uideal[3].imag = 0;

	return fidelity(dim, U, Uideal);
}
