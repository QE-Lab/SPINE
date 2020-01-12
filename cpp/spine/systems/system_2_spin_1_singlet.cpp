/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "system_2_spin_1_singlet.h"

using namespace spine;

spine::systems::system_2_spin_1_singlet::system_2_spin_1_singlet()
	: spine::systems::system_spin(2)
{
	
}

spine::systems::system_2_spin_1_singlet::~system_2_spin_1_singlet()
{

}

int spine::systems::system_2_spin_1_singlet::getIndex(unsigned int * states)
{
	if ((states[0] == STATE_0) && (states[1] == STATE_0))
		return 0;
	if ((states[0] == STATE_0) && (states[1] == STATE_1))
		return 1;
	if ((states[0] == STATE_1) && (states[1] == STATE_0))
		return 2;
	if ((states[0] == STATE_1) && (states[1] == STATE_1))
		return 3;
	if ((states[0] == STATE_N) && (states[1] == STATE_S))
		return 4;
	return -1;
}

vector<unsigned int> * spine::systems::system_2_spin_1_singlet::getIndexMeasurement(unsigned int dot, unsigned int state)
{
	if (dot == 0) {
		if (state == STATE_0)
			return &measurement_0_0;
		if (state == STATE_1)
			return &measurement_0_1;
		if (state == STATE_N)
			return &measurement_0_N;
	}
	if (dot == 1) {
		if (state == STATE_0)
			return &measurement_1_0;
		if (state == STATE_1)
			return &measurement_1_1;
		if (state == STATE_S)
			return &measurement_1_S;
	}
	return nullptr;
}

unsigned int spine::systems::system_2_spin_1_singlet::getDimension()
{
	return 5;
}

void spine::systems::system_2_spin_1_singlet::updateHamiltonian(realnum * H)
{
	// Single microwave driveline
	realnum microwave_control = this->microwave_control[0] + this->microwave_control[1];

	// Common Rabi frequency in simple Hamiltonian
	realnum rabi_frequency = this->rabi_frequency[0];

	// All tunnel rates are the same for a double dot
	realnum tunnel_control = this->tunnel_control[0];

	// Only singlet for dot 2
	realnum charging_energy = this->charging_energy[1];
	realnum detuning_control = this->detuning_control[1];

	H[0] = -(larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
	H[1] = microwave_control * rabi_frequency;
	H[2] = microwave_control * rabi_frequency;
	H[3] = 0;
	H[4] = 0;
	H[5] = microwave_control * rabi_frequency;
	H[6] = -(larmor_frequency[0] - larmor_frequency[1]) / ((realnum) 2.0);
	H[7] = 0;
	H[8] = microwave_control * rabi_frequency;
	H[9] = tunnel_control;
	H[10] = microwave_control * rabi_frequency;
	H[11] = 0;
	H[12] = +(larmor_frequency[0] - larmor_frequency[1]) / ((realnum) 2.0);
	H[13] = microwave_control * rabi_frequency;
	H[14] = -tunnel_control;
	H[15] = 0;
	H[16] = microwave_control * rabi_frequency;
	H[17] = microwave_control * rabi_frequency;
	H[18] = +(larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
	H[19] = 0;
	H[20] = 0;
	H[21] = tunnel_control;
	H[22] = -tunnel_control;
	H[23] = 0;
	H[24] = charging_energy - detuning_control;
}
