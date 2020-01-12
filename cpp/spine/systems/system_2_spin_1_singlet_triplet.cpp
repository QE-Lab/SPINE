/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "system_2_spin_1_singlet_triplet.h"

using namespace spine;

spine::systems::system_2_spin_1_singlet_triplet::system_2_spin_1_singlet_triplet()
	: spine::systems::system_spin(2)
{
	
}

spine::systems::system_2_spin_1_singlet_triplet::~system_2_spin_1_singlet_triplet()
{

}

int spine::systems::system_2_spin_1_singlet_triplet::getIndex(unsigned int * states)
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
	if ((states[0] == STATE_N) && (states[1] == STATE_T0))
		return 5;
	if ((states[0] == STATE_N) && (states[1] == STATE_TP))
		return 6;
	if ((states[0] == STATE_N) && (states[1] == STATE_TM))
		return 7;
	return -1;
}

vector<unsigned int> * spine::systems::system_2_spin_1_singlet_triplet::getIndexMeasurement(unsigned int dot, unsigned int state)
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
		if (state == STATE_T0)
			return &measurement_1_T0;
		if (state == STATE_TP)
			return &measurement_1_TP;
		if (state == STATE_TM)
			return &measurement_1_TM;
	}
	return nullptr;
}

unsigned int spine::systems::system_2_spin_1_singlet_triplet::getDimension()
{
	return 8;
}

void spine::systems::system_2_spin_1_singlet_triplet::updateHamiltonian(realnum * H)
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
	realnum singlet_triplet_energy = this->singlet_triplet_energy[1];

	H[0] = -(larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
	H[1] = microwave_control * rabi_frequency;
	H[2] = microwave_control * rabi_frequency;
	H[3] = 0;
	H[4] = sqrt(2) * tunnel_control;
	H[5] = 0;
	H[6] = 0;
	H[7] = 0;
	H[8] = microwave_control*rabi_frequency;
	H[9] = -(larmor_frequency[0] - larmor_frequency[1]) / ((realnum) 2.0);
	H[10] = 0;
	H[11] = microwave_control * rabi_frequency;
	H[12] = 0;
	H[13] = sqrt(2) * tunnel_control;
	H[14] = 0;
	H[15] = 0;
	H[16] = microwave_control*rabi_frequency;
	H[17] = 0;
	H[18] = +(larmor_frequency[0] - larmor_frequency[1]) / ((realnum) 2.0);
	H[19] = microwave_control * rabi_frequency;
	H[20] = 0;
	H[21] = 0;
	H[22] = sqrt(2) * tunnel_control;
	H[23] = 0;
	H[24] = 0;
	H[25] = microwave_control * rabi_frequency;
	H[26] = microwave_control * rabi_frequency;
	H[27] = +(larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
	H[28] = 0;
	H[29] = 0;
	H[30] = 0;
	H[31] = sqrt(2) * tunnel_control;
	H[32] = sqrt(2) * tunnel_control;
	H[33] = 0;
	H[34] = 0;
	H[35] = 0;
	H[36] = charging_energy - detuning_control + singlet_triplet_energy - (larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
	H[37] = 0;
	H[38] = 0;
	H[39] = 0;
	H[40] = 0;
	H[41] = sqrt(2) * tunnel_control;
	H[42] = 0;
	H[43] = 0;
	H[44] = 0;
	H[45] = charging_energy - detuning_control + singlet_triplet_energy / ((realnum) 2.0);
	H[46] = singlet_triplet_energy / 2.0;
	H[47] = 0;
	H[48] = 0;
	H[49] = 0;
	H[50] = sqrt(2) * tunnel_control;
	H[51] = 0;
	H[52] = 0;
	H[53] = singlet_triplet_energy / 2.0;
	H[54] = charging_energy - detuning_control + singlet_triplet_energy / ((realnum) 2.0);
	H[55] = 0;
	H[56] = 0;
	H[57] = 0;
	H[58] = 0;
	H[59] = sqrt(2) * tunnel_control;
	H[60] = 0;
	H[61] = 0;
	H[62] = 0;
	H[63] = charging_energy - detuning_control + singlet_triplet_energy + (larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
}
