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
		class system_dispersive_readout {

		protected:

			// Boltzmann constant(eV)
			const realnum kB = (realnum) 8.62e-5;

			// Hamiltonian Properties
			realnum electron_temperature;
			realnum tunnel_rate;
			realnum lever_arm;

			// Hamiltonian Control
			realnum gate_voltage;

			// Plot handles
#ifdef PLOT
			unsigned int points = 100;
			double xmin = 0;
			double xmax = 0;
			Plot * p = nullptr;
#endif

		public:
			system_dispersive_readout();
			~system_dispersive_readout();

			unsigned int getDimension();
			void initialize(complex * state, realnum P);
			realnum measure(complex * state);
			void updateHamiltonian(complex * H);

#ifdef PLOT
			void plotSetup(unsigned int points = 100, double xmin = 0, double xmax = 0);
			void plotAddU(complex * U, realnum t = 0, bool plot_lab = false);
			void plotAdd(complex * state, realnum t = 0, bool plot_lab = false);
			void plot();
#endif

			void setElectronTemperature(realnum value);
			realnum getElectronTemperature();
			void setTunnelRate(realnum value);
			realnum getTunnelRate();
			void setLeverArm(realnum value);
			realnum getLeverArm();

			void setGateVoltage(realnum value);
		};
	}
}
