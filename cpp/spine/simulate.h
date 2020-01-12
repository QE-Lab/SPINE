/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#pragma once

#include "math.h"

using namespace spine::math;

namespace spine
{
	// Function prototypes
	void simulate(unsigned int dim,
		bool(*inHamiltonian)(realnum * H, realnum * timestep),
		void(*outOperation)(complex * U),
		void(*solver)(unsigned int dim, realnum * H, complex * dU));
	void simulate(unsigned int dim,
		bool(*inHamiltonian)(realnum * H, realnum * timestep),
		void(*outOperation)(complex * U),
		void(*solver)(unsigned int dim, realnum * H, complex * dU),
		complex * state);
	void simulate(unsigned int dim,
		bool(*inHamiltonian)(complex * H, realnum * timestep),
		void(*outOperation)(complex * U),
		void(*solver)(unsigned int dim, complex * H, complex * dU));
	void simulate(unsigned int dim,
		bool(*inHamiltonian)(complex * H, realnum * timestep),
		void(*outOperation)(complex * U),
		void(*solver)(unsigned int dim, complex * H, complex * dU),
		complex * state);
}
