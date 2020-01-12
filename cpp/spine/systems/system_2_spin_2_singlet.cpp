/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "system_2_spin_2_singlet.h"

using namespace spine;

spine::systems::system_2_spin_2_singlet::system_2_spin_2_singlet()
	: spine::systems::system_spin(2)
{
	
}

spine::systems::system_2_spin_2_singlet::~system_2_spin_2_singlet()
{

}

int spine::systems::system_2_spin_2_singlet::getIndex(unsigned int * states)
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
	if ((states[0] == STATE_S) && (states[1] == STATE_N))
		return 5;
	return -1;
}

vector<unsigned int> * spine::systems::system_2_spin_2_singlet::getIndexMeasurement(unsigned int dot, unsigned int state)
{
	if (dot == 0) {
		if (state == STATE_0)
			return &measurement_0_0;
		if (state == STATE_1)
			return &measurement_0_1;
		if (state == STATE_N)
			return &measurement_0_N;
		if (state == STATE_S)
			return &measurement_0_S;
	}
	if (dot == 1) {
		if (state == STATE_0)
			return &measurement_1_0;
		if (state == STATE_1)
			return &measurement_1_1;
		if (state == STATE_S)
			return &measurement_1_S;
		if (state == STATE_N)
			return &measurement_1_N;
	}
	return nullptr;
}

unsigned int spine::systems::system_2_spin_2_singlet::getDimension()
{
	return 6;
}

void spine::systems::system_2_spin_2_singlet::updateHamiltonian(realnum * H)
{
	// Single microwave driveline
	realnum microwave_control = this->microwave_control[0] + this->microwave_control[1];

	// Common Rabi frequency and charging energy in simple Hamiltonian
	realnum rabi_frequency = this->rabi_frequency[0];
	realnum charging_energy = this->charging_energy[0];

	// All tunnel rates are the same for a double dot
	realnum tunnel_control = this->tunnel_control[0];

	// Relative detuning
	realnum detuning_control = this->detuning_control[1] - this->detuning_control[0];

	H[0] = -(larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
	H[1] = microwave_control * rabi_frequency;
	H[2] = microwave_control * rabi_frequency;
	H[3] = 0;
	H[4] = 0;
	H[5] = 0;
	H[6] = microwave_control * rabi_frequency;
	H[7] = -(larmor_frequency[0] - larmor_frequency[1]) / ((realnum) 2.0);
	H[8] = 0;
	H[9] = microwave_control * rabi_frequency;
	H[10] = tunnel_control;
	H[11] = tunnel_control;
	H[12] = microwave_control * rabi_frequency;
	H[13] = 0;
	H[14] = +(larmor_frequency[0] - larmor_frequency[1]) / ((realnum) 2.0);
	H[15] = microwave_control * rabi_frequency;
	H[16] = -tunnel_control;
	H[17] = -tunnel_control;
	H[18] = 0;
	H[19] = microwave_control * rabi_frequency;
	H[20] = microwave_control * rabi_frequency;
	H[21] = +(larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
	H[22] = 0;
	H[23] = 0;
	H[24] = 0;
	H[25] = tunnel_control;
	H[26] = -tunnel_control;
	H[27] = 0;
	H[28] = charging_energy - detuning_control;
	H[29] = 0;
	H[30] = 0;
	H[31] = tunnel_control;
	H[32] = -tunnel_control;
	H[33] = 0;
	H[34] = 0;
	H[35] = charging_energy + detuning_control;
}
