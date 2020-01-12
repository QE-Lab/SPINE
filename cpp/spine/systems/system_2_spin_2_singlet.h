/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include "../math.h"
#include "system_spin.h"

using namespace std;
using namespace spine::math;

namespace spine
{
	namespace systems
	{
		class system_2_spin_2_singlet : public spine::systems::system_spin {

		private:
			vector<unsigned int> measurement_0_0{ 0, 1 };
			vector<unsigned int> measurement_0_1{ 2, 3 };
			vector<unsigned int> measurement_0_N{ 4 };
			vector<unsigned int> measurement_0_S{ 5 };
			vector<unsigned int> measurement_1_0{ 0, 2 };
			vector<unsigned int> measurement_1_1{ 1, 3 };
			vector<unsigned int> measurement_1_S{ 4 };
			vector<unsigned int> measurement_1_N{ 5 };

		public:
			system_2_spin_2_singlet();
			~system_2_spin_2_singlet();

			int getIndex(unsigned int * states);
			vector<unsigned int> * getIndexMeasurement(unsigned int dot, unsigned int state);
			unsigned int getDimension();
			void updateHamiltonian(realnum * H);
		};
	}
}
