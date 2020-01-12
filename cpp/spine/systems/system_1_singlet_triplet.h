/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include "../math.h"
#ifdef PLOT
#include "../plot.h"
#endif

using namespace std;
using namespace spine::math;

namespace spine
{
	namespace systems
	{
		class system_1_singlet_triplet {

		protected:

			// Hamiltonian Properties
			realnum magnetic_gradient;

			// Hamiltonian Control
			realnum exchange_interaction;

			// Plot handles
#ifdef PLOT
			unsigned int points = 100;
			double xmin = 0;
			double xmax = 0;
			Plot * p = nullptr;
#endif

		public:
			system_1_singlet_triplet();
			~system_1_singlet_triplet();

			unsigned int getDimension();
			void initialize(complex * state);
			void measure(complex * state, realnum * Px, realnum * Py, realnum * Pz);
			void updateHamiltonian(realnum * H);

#ifdef PLOT
			void plotSetup(unsigned int points = 100, double xmin = 0, double xmax = 0);
			void plotAddU(complex * U, realnum t = 0, bool plot_lab = false);
			void plotAdd(complex * state, realnum t = 0, bool plot_lab = false);
			void plot();
#endif

			void setMagneticGradient(realnum value);
			realnum getMagneticGradient();

			void setExchangeInteraction(realnum value);
		};
	}
}
