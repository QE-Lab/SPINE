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
		class system_n_spin_n_singlet : public spine::systems::system_spin {

		private:

			typedef struct {
				unsigned int x;
				unsigned int y;
			} storage_location;

			// Information about the basis
			unsigned int size;
			unsigned int * storage_lookup;
			unsigned int storage_size;
			vector<unsigned int> * probability_i;

			// Information for the Hamiltonian built-up
			vector<unsigned int> * charging_energies_i;
			vector<unsigned int> * larmor_frequencies_ip;
			vector<unsigned int> * larmor_frequencies_im;
			vector<storage_location> * tunnel_couplings_ip;
			vector<storage_location> * tunnel_couplings_im;
			vector<storage_location> * microwave_control_i;
			vector<unsigned int> * detuning_control_ip;
			vector<unsigned int> * detuning_control_im;

		private:
			bool isValid(unsigned int index);

		public:
			system_n_spin_n_singlet(unsigned int dots);
			~system_n_spin_n_singlet();

			int getIndex(unsigned int * states);
			vector<unsigned int> * getIndexMeasurement(unsigned int dot, unsigned int state);
			unsigned int getDimension();
			void updateHamiltonian(realnum * H);
		};
	}
}
