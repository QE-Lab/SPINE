/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "system_1_spin_rwa.h"

using namespace spine;

spine::systems::system_1_spin_rwa::system_1_spin_rwa()
	: spine::systems::system_spin(1)
{

}

spine::systems::system_1_spin_rwa::~system_1_spin_rwa()
{

}

int spine::systems::system_1_spin_rwa::getIndex(unsigned int * states)
{
	if (states[0] == STATE_0)
		return 0;
	if (states[0] == STATE_1)
		return 1;
	return -1;
}

vector<unsigned int> * spine::systems::system_1_spin_rwa::getIndexMeasurement(unsigned int dot, unsigned int state)
{
	if (dot == 0) {
		if (state == STATE_0)
			return &measurement_0_0;
		if (state == STATE_1)
			return &measurement_0_1;
	}
	return nullptr;
}

unsigned int spine::systems::system_1_spin_rwa::getDimension()
{
	return 2;
}

void spine::systems::system_1_spin_rwa::updateHamiltonian(complex * H)
{
	H[0].real = (microwave_frequency - larmor_frequency[0]) / ((realnum) 2.0);
	H[0].imag = 0;
	H[1].real = (rabi_frequency[0] * microwave_amplitude) / ((realnum) 2.0) * (realnum) cos(microwave_phase);
	H[1].imag = (rabi_frequency[0] * microwave_amplitude) / ((realnum) 2.0) * (realnum) sin(microwave_phase);
	H[2].real = H[1].real;
	H[2].imag = -H[1].imag;
	H[3].real = -H[0].real;
	H[3].imag = 0;
}

void spine::systems::system_1_spin_rwa::setMicrowaveControl(realnum value)
{
	assert(false);
}

void spine::systems::system_1_spin_rwa::setMicrowaveControl(unsigned int dot, realnum value)
{
	assert(false);
}

void spine::systems::system_1_spin_rwa::setMicrowaveFrequency(realnum value)
{
	microwave_frequency = value;
}

realnum spine::systems::system_1_spin_rwa::getMicrowaveFrequency()
{
	return microwave_frequency;
}

void spine::systems::system_1_spin_rwa::setMicrowaveAmplitude(realnum value)
{
	microwave_amplitude = value;
}

void spine::systems::system_1_spin_rwa::setMicrowavePhase(realnum value)
{
	microwave_phase = value;
}
