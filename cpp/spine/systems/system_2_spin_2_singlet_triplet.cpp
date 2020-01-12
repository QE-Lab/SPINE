/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "system_2_spin_2_singlet_triplet.h"

using namespace spine;

spine::systems::system_2_spin_2_singlet_triplet::system_2_spin_2_singlet_triplet()
	: spine::systems::system_spin(2)
{
	
}

spine::systems::system_2_spin_2_singlet_triplet::~system_2_spin_2_singlet_triplet()
{

}

int spine::systems::system_2_spin_2_singlet_triplet::getIndex(unsigned int * states)
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
	if ((states[0] == STATE_S) && (states[1] == STATE_N))
		return 8;
	if ((states[0] == STATE_T0) && (states[1] == STATE_N))
		return 9;
	if ((states[0] == STATE_TP) && (states[1] == STATE_N))
		return 10;
	if ((states[0] == STATE_TM) && (states[1] == STATE_N))
		return 11;
	return -1;
}

vector<unsigned int> * spine::systems::system_2_spin_2_singlet_triplet::getIndexMeasurement(unsigned int dot, unsigned int state)
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
		if (state == STATE_T0)
			return &measurement_0_T0;
		if (state == STATE_TP)
			return &measurement_0_TP;
		if (state == STATE_TM)
			return &measurement_0_TM;
	}
	if (dot == 1) {
		if (state == STATE_0)
			return &measurement_1_0;
		if (state == STATE_1)
			return &measurement_1_1;
		if (state == STATE_N)
			return &measurement_1_N;
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

unsigned int spine::systems::system_2_spin_2_singlet_triplet::getDimension()
{
	return 12;
}

void spine::systems::system_2_spin_2_singlet_triplet::updateHamiltonian(realnum * H)
{
	// Single microwave driveline
	realnum microwave_control = this->microwave_control[0] + this->microwave_control[1];

	// Common Rabi frequency, charging energy and singlet-triplet splitting in simple Hamiltonian
	realnum rabi_frequency = this->rabi_frequency[0];
	realnum charging_energy = this->charging_energy[0];
	realnum singlet_triplet_energy = this->singlet_triplet_energy[0];

	// All tunnel rates are the same for a double dot
	realnum tunnel_control = this->tunnel_control[0];

	// Relative detuning
	realnum detuning_control = this->detuning_control[1] - this->detuning_control[0];
	
	H[0] = -(larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
	H[1] = microwave_control * rabi_frequency;
	H[2] = microwave_control * rabi_frequency;
	H[3] = 0;
	H[4] = sqrt(2) * tunnel_control;
	H[5] = 0;
	H[6] = 0;
	H[7] = 0;
	H[8] = sqrt(2) * tunnel_control;
	H[9] = 0;
	H[10] = 0;
	H[11] = 0;
	H[12] = microwave_control * rabi_frequency;
	H[13] = -(larmor_frequency[0] - larmor_frequency[1]) / ((realnum) 2.0);
	H[14] = 0;
	H[15] = microwave_control * rabi_frequency;
	H[16] = 0;
	H[17] = sqrt(2) * tunnel_control;
	H[18] = 0;
	H[19] = 0;
	H[20] = 0;
	H[21] = sqrt(2) * tunnel_control;
	H[22] = 0;
	H[23] = 0;
	H[24] = microwave_control * rabi_frequency;
	H[25] = 0;
	H[26] = +(larmor_frequency[0] - larmor_frequency[1]) / ((realnum) 2.0);
	H[27] = microwave_control * rabi_frequency;
	H[28] = 0;
	H[29] = 0;
	H[30] = sqrt(2) * tunnel_control;
	H[31] = 0;
	H[32] = 0;
	H[33] = 0;
	H[34] = sqrt(2) * tunnel_control;
	H[35] = 0;
	H[36] = 0;
	H[37] = microwave_control * rabi_frequency;
	H[38] = microwave_control * rabi_frequency;
	H[39] = +(larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
	H[40] = 0;
	H[41] = 0;
	H[42] = 0;
	H[43] = sqrt(2) * tunnel_control;
	H[44] = 0;
	H[45] = 0;
	H[46] = 0;
	H[47] = sqrt(2) * tunnel_control;
	H[48] = sqrt(2) * tunnel_control;
	H[49] = 0;
	H[50] = 0;
	H[51] = 0;
	H[52] = charging_energy - detuning_control + singlet_triplet_energy - (larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
	H[53] = 0;
	H[54] = 0;
	H[55] = 0;
	H[56] = 0;
	H[57] = 0;
	H[58] = 0;
	H[59] = 0;
	H[60] = 0;
	H[61] = sqrt(2) * tunnel_control;
	H[62] = 0;
	H[63] = 0;
	H[64] = 0;
	H[65] = charging_energy - detuning_control + singlet_triplet_energy / ((realnum) 2.0);
	H[66] = singlet_triplet_energy / ((realnum) 2.0);
	H[67] = 0;
	H[68] = 0;
	H[69] = 0;
	H[70] = 0;
	H[71] = 0;
	H[72] = 0;
	H[73] = 0;
	H[74] = sqrt(2) * tunnel_control;
	H[75] = 0;
	H[76] = 0;
	H[77] = singlet_triplet_energy / ((realnum) 2.0);
	H[78] = charging_energy - detuning_control + singlet_triplet_energy / ((realnum) 2.0);
	H[79] = 0;
	H[80] = 0;
	H[81] = 0;
	H[82] = 0;
	H[83] = 0;
	H[84] = 0;
	H[85] = 0;
	H[86] = 0;
	H[87] = sqrt(2) * tunnel_control;
	H[88] = 0;
	H[89] = 0;
	H[90] = 0;
	H[91] = charging_energy - detuning_control + singlet_triplet_energy + (larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
	H[92] = 0;
	H[93] = 0;
	H[94] = 0;
	H[95] = 0;
	H[96] = sqrt(2) * tunnel_control;
	H[97] = 0;
	H[98] = 0;
	H[99] = 0;
	H[100] = 0;
	H[101] = 0;
	H[102] = 0;
	H[103] = 0;
	H[104] = charging_energy + detuning_control + singlet_triplet_energy - (larmor_frequency[0] + larmor_frequency[1]) / ((realnum) 2.0);
	H[105] = 0;
	H[106] = 0;
	H[107] = 0;
	H[108] = 0;
	H[109] = sqrt(2) * tunnel_control;
	H[110] = 0;
	H[111] = 0;
	H[112] = 0;
	H[113] = 0;
	H[114] = 0;
	H[115] = 0;
	H[116] = 0;
	H[117] = charging_energy + detuning_control + singlet_triplet_energy / ((realnum) 2.0);
	H[118] = singlet_triplet_energy / ((realnum) 2.0);
	H[119] = 0;
	H[120] = 0;
	H[121] = 0;
	H[122] = sqrt(2) * tunnel_control;
	H[123] = 0;
	H[124] = 0;
	H[125] = 0;
	H[126] = 0;
	H[127] = 0;
	H[128] = 0;
	H[129] = singlet_triplet_energy / ((realnum) 2.0);
	H[130] = charging_energy + detuning_control + singlet_triplet_energy / ((realnum) 2.0);
	H[131] = 0;
	H[132] = 0;
	H[133] = 0;
	H[134] = 0;
	H[135] = sqrt(2) * tunnel_control;
	H[136] = 0;
	H[137] = 0;
	H[138] = 0;
	H[139] = 0;
	H[140] = 0;
	H[141] = 0;
	H[142] = 0;
	H[143] = charging_energy + detuning_control + singlet_triplet_energy + (larmor_frequency[0] + larmor_frequency[1]) / 2.0;
}
