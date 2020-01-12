/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "system_1_spin.h"

using namespace spine;

spine::systems::system_1_spin::system_1_spin()
	: spine::systems::system_spin(1)
{
	
}

spine::systems::system_1_spin::~system_1_spin()
{

}

int spine::systems::system_1_spin::getIndex(unsigned int * states)
{
	if (states[0] == STATE_0)
		return 0;
	if (states[0] == STATE_1)
		return 1;
	return -1;
}

vector<unsigned int> * spine::systems::system_1_spin::getIndexMeasurement(unsigned int dot, unsigned int state)
{
	if (dot == 0) {
		if (state == STATE_0)
			return &measurement_0_0;
		if (state == STATE_1)
			return &measurement_0_1;
	}
	return nullptr;
}

unsigned int spine::systems::system_1_spin::getDimension()
{
	return 2;
}

void spine::systems::system_1_spin::updateHamiltonian(realnum * H)
{
	H[0] = -larmor_frequency[0] / ((realnum) 2.0);
	H[1] = rabi_frequency[0] * microwave_control[0];
	H[2] = H[1];
	H[3] = -H[0];
}
