/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#pragma once

#include <cassert>
#include "../math.h"

using namespace spine::math;

namespace spine
{
	namespace solvers
	{
		void solver_analytical_xz(unsigned int dim, realnum * dtH, complex * dU);
	}
}
