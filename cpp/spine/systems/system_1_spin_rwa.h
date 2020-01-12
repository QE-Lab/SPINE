/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <cassert>
#include <vector>
#include "../math.h"
#include "system_spin.h"

using namespace std;
using namespace spine::math;

namespace spine
{
	namespace systems
	{
		class system_1_spin_rwa : public spine::systems::system_spin {

		protected:

			// Hamiltonian Properties
			realnum microwave_frequency;

			// Hamiltonian Control
			realnum microwave_amplitude;
			realnum microwave_phase;

		private:
			vector<unsigned int> measurement_0_0{ 0 };
			vector<unsigned int> measurement_0_1{ 1 };

		public:
			system_1_spin_rwa();
			~system_1_spin_rwa();

			int getIndex(unsigned int * states);
			vector<unsigned int> * getIndexMeasurement(unsigned int dot, unsigned int state);
			unsigned int getDimension();
			void updateHamiltonian(complex * H);

			void setMicrowaveControl(realnum value);
			void setMicrowaveControl(unsigned int dot, realnum value);

			void setMicrowaveFrequency(realnum value);
			realnum getMicrowaveFrequency();

			void setMicrowaveAmplitude(realnum value);
			void setMicrowavePhase(realnum value);
		};
	}
}
