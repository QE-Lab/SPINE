/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <cwchar>
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
		class system_spin {

		protected:

			// Number of dots
			unsigned int dots;

			// Hamiltonian Properties
			realnum * larmor_frequency;
			realnum * rabi_frequency;
			realnum * charging_energy;
			realnum * singlet_triplet_energy;

			// Hamiltonian Control
			realnum * microwave_control;
			realnum * detuning_control;
			realnum * tunnel_control;

			// Plot handles
#ifdef PLOT
			unsigned int points = 100;
			double xmin = 0;
			double xmax = 0;
			Plot ** p = nullptr;
#endif

		public:

			// Supported spin states :
			const unsigned int STATE_N = 0;			// 0 electrons
			const unsigned int STATE_0 = 1;			// 1 electron in ground state
			const unsigned int STATE_1 = 2;			// 1 electron in excited state
			const unsigned int STATE_S = 3;			// 2 electrons in singlet configuration(or in 00 spin state, depending on the Hamiltonian)
			const unsigned int STATE_T0 = 4;		// 2 electrons in triplet(0) configuration(or in 01 spin state, depending on the Hamiltonian)
			const unsigned int STATE_TP = 5;		// 2 electrons in triplet(+) configuration(or in 10 spin state, depending on the Hamiltonian)
			const unsigned int STATE_TM = 6;		// 2 electrons in triplet(-) configuration(or in 11 spin state, depending on the Hamiltonian)

		public:

			virtual int getIndex(unsigned int * states) = 0;
			virtual vector<unsigned int> * getIndexMeasurement(unsigned int dot, unsigned int state) = 0;
			virtual unsigned int getDimension() = 0;

		public:
			system_spin(unsigned int dots);
			~system_spin();

			int initializeIndex();
			void initialize(complex * state);
			void measure(complex * state, realnum t, realnum * Px, realnum * Py, realnum * Pz);
			void measure(complex * state, realnum * Px, realnum * Py, realnum * Pz);
			void measureST(complex * state, realnum * Pn, realnum * Ps, realnum * Pt);

#ifdef PLOT
			void plotSetup(unsigned int points = 100, double xmin = 0, double xmax = 0);
			void plotAddU(complex * U, realnum t = 0, bool plot_lab = false);
			void plotAdd(complex * state, realnum t = 0, bool plot_lab = false);
			void plot();
#endif

			void setLarmorFrequency(realnum value);
			void setLarmorFrequency(unsigned int dot, realnum value);
			realnum getLarmorFrequency();
			realnum getLarmorFrequency(unsigned int dot);

			void setRabiFrequency(realnum value);
			void setRabiFrequency(unsigned int dot, realnum value);
			realnum getRabiFrequency();
			realnum getRabiFrequency(unsigned int dot);

			void setChargingEnergy(realnum value);
			void setChargingEnergy(unsigned int dot, realnum value);
			realnum getChargingEnergy();
			realnum getChargingEnergy(unsigned int dot);

			void setSingletTripletEnergy(realnum value);
			void setSingletTripletEnergy(unsigned int dot, realnum value);
			realnum getSingletTripletEnergy();
			realnum getSingletTripletEnergy(unsigned int dot);

			void setMicrowaveControl(realnum value);
			void setMicrowaveControl(unsigned int dot, realnum value);

			void setDetuningControl(realnum value);
			void setDetuningControl(unsigned int dot, realnum value);

			void setTunnelControl(realnum value);
			void setTunnelControl(unsigned int dota, unsigned int dotb, realnum value);
		};
	}
}
