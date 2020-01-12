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
		class system_2_spin_2_singlet_triplet : public spine::systems::system_spin {

		private:
			vector<unsigned int> measurement_0_0{ 0, 1 };
			vector<unsigned int> measurement_0_1{ 2, 3 };
			vector<unsigned int> measurement_0_N{ 4, 5, 6, 7 };
			vector<unsigned int> measurement_0_S{ 8 };
			vector<unsigned int> measurement_0_T0{ 9 };
			vector<unsigned int> measurement_0_TP{ 10 };
			vector<unsigned int> measurement_0_TM{ 11 };
			vector<unsigned int> measurement_1_0{ 0, 2 };
			vector<unsigned int> measurement_1_1{ 1, 3 };
			vector<unsigned int> measurement_1_N{ 8, 9, 10, 11 };
			vector<unsigned int> measurement_1_S{ 4 };
			vector<unsigned int> measurement_1_T0{ 5 };
			vector<unsigned int> measurement_1_TP{ 6 };
			vector<unsigned int> measurement_1_TM{ 7 };

		public:
			system_2_spin_2_singlet_triplet();
			~system_2_spin_2_singlet_triplet();

			int getIndex(unsigned int * states);
			vector<unsigned int> * getIndexMeasurement(unsigned int dot, unsigned int state);
			unsigned int getDimension();
			void updateHamiltonian(realnum * H);
		};
	}
}
