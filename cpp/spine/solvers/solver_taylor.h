/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#pragma once

#include "../math.h"

using namespace spine::math;

namespace spine
{
	namespace solvers
	{
		void solver_taylor(unsigned int dim, realnum * dtH, complex * dU);
		void solver_taylor(unsigned int dim, complex * dtH, complex * dU);
	}
}
