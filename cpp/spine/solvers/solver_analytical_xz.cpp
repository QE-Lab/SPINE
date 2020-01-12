/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "solver_analytical_xz.h"

using namespace spine;

void spine::solvers::solver_analytical_xz(unsigned int dim, realnum * dtH, complex * dU)
{
	assert(dim == 2);

	// Note: analytically solved, assumes the following form for dtH
	// dtH = [a,  b
	//        b, -a]
	realnum a = dtH[0];
	realnum b = dtH[1];

	// dU = expm(-1i * dtH)
#ifdef DOUBLE_PRECISION

	realnum tmp = sqrt(a*a + b * b);
	dU[0].real = cos(tmp);
	dU[0].imag = -a / tmp * sin(tmp);
	dU[1].imag = -b / tmp * sin(tmp);

#else

	realnum tmp = sqrtf(a*a + b * b);
	dU[0].real = cosf(tmp);
	dU[0].imag = -a / tmp * sinf(tmp);
	dU[1].imag = -b / tmp * sinf(tmp);

#endif

	dU[1].real = 0;
	dU[2].real = 0;
	dU[2].imag = dU[1].imag;
	dU[3].real = dU[0].real;
	dU[3].imag = -dU[0].imag;
}
