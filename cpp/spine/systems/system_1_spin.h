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
		class system_1_spin : public spine::systems::system_spin {

		private:
			vector<unsigned int> measurement_0_0{ 0 };
			vector<unsigned int> measurement_0_1{ 1 };

		public:
			system_1_spin();
			~system_1_spin();

			int getIndex(unsigned int * states);
			vector<unsigned int> * getIndexMeasurement(unsigned int dot, unsigned int state);
			unsigned int getDimension();
			void updateHamiltonian(realnum * H);
		};
	}
}
