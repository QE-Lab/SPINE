/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "system_dispersive_readout.h"

using namespace spine;

spine::systems::system_dispersive_readout::system_dispersive_readout()
{
	
}

spine::systems::system_dispersive_readout::~system_dispersive_readout()
{
#ifdef PLOT
	delete p;
#endif
}

unsigned int spine::systems::system_dispersive_readout::getDimension()
{
	return 2;
}

void spine::systems::system_dispersive_readout::initialize(complex * state, realnum P)
{
	state[0].real = 1;
	state[0].imag = 0;
	state[1].real = P;
	state[1].imag = 0;
}

realnum spine::systems::system_dispersive_readout::measure(complex * state)
{
	return state[1].real;
}

void spine::systems::system_dispersive_readout::updateHamiltonian(complex * H)
{
	H[0] = { 0, 0 };
	H[1] = { 0, 0 };
	H[2].real = 0;
	H[2].imag = tunnel_rate / (1 + exp(lever_arm * gate_voltage / kB / electron_temperature));
	H[3].real = 0;
	H[3].imag = -tunnel_rate;
}

#ifdef PLOT

void spine::systems::system_dispersive_readout::plotSetup(unsigned int points, double xmin, double xmax)
{
	this->points = points;
	this->xmin = xmin;
	this->xmax = xmax;
}

void spine::systems::system_dispersive_readout::plotAddU(complex * U, realnum t, bool plot_lab)
{
	unsigned int dim = getDimension();
	complex * statetmp = allocComplex(dim);
	complex * state = allocComplex(dim);

	// Update the overall operation
	initialize(state, 0.0);
	initCopy(dim, statetmp, state);
	multiplyVector(dim, U, statetmp, state);
	plotAdd(state, t, plot_lab);

	freeComplex(state);
	freeComplex(statetmp);
}

void spine::systems::system_dispersive_readout::plotAdd(complex * state, realnum t, bool plot_lab)
{
	realnum P;
	P = measure(state);

	if (p == nullptr)
	{
		p = new Plot(L"Probability", points);
		p->add(0, 0, 255, xmin, xmax, 0, 1);
	}
	p->get(0)->add(t, P);
}

void spine::systems::system_dispersive_readout::plot()
{
	if (p != nullptr)
		p->redraw();
}

#endif

void spine::systems::system_dispersive_readout::setElectronTemperature(realnum value)
{
	electron_temperature = value;
}

realnum spine::systems::system_dispersive_readout::getElectronTemperature()
{
	return electron_temperature;
}

void spine::systems::system_dispersive_readout::setTunnelRate(realnum value)
{
	tunnel_rate = value;
}

realnum spine::systems::system_dispersive_readout::getTunnelRate()
{
	return tunnel_rate;
}
void spine::systems::system_dispersive_readout::setLeverArm(realnum value)
{
	lever_arm = value;
}

realnum spine::systems::system_dispersive_readout::getLeverArm()
{
	return lever_arm;
}
void spine::systems::system_dispersive_readout::setGateVoltage(realnum value)
{
	gate_voltage = value;
}
