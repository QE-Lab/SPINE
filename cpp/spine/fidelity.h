/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#pragma once

#include <cassert>
#include "math.h"

using namespace spine::math;

namespace spine
{
	// Function prototypes
	realnum fidelity(unsigned int dim, complex * U, complex * Uideal);
	realnum fidelity(unsigned int dim, complex * U, realnum theta, realnum phi);
}
