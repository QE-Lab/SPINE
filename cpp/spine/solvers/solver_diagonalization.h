/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#pragma once

#include "../math.h"

using namespace spine::math;

#ifdef MKL

namespace spine
{
	namespace solvers
	{
		void solver_diagonalization(unsigned int dim, realnum * dtH, complex * dU);
		void solver_diagonalization(unsigned int dim, complex * dtH, complex * dU);
	}
}

#endif
