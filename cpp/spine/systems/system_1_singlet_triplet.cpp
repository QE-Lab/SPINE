/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "system_1_singlet_triplet.h"

using namespace spine;

spine::systems::system_1_singlet_triplet::system_1_singlet_triplet()
{
	magnetic_gradient = 0;
	exchange_interaction = 0;
}

spine::systems::system_1_singlet_triplet::~system_1_singlet_triplet()
{
#ifdef PLOT
	delete p;
#endif
}

unsigned int spine::systems::system_1_singlet_triplet::getDimension()
{
	return 2;
}

void spine::systems::system_1_singlet_triplet::initialize(complex * state)
{
	state[0].real = 1;
	state[0].imag = 0;
	state[1].real = 0;
	state[1].imag = 0;
}

void spine::systems::system_1_singlet_triplet::measure(complex * state, realnum * Px, realnum * Py, realnum * Pz)
{
	complex tmp;

	tmp.real = state[0].real + state[1].real;
	tmp.imag = state[0].imag + state[1].imag;
	Px[0] = abs2(tmp) / ((realnum) 2.0);

	tmp.real = state[0].real - state[1].imag;
	tmp.imag = state[0].imag + state[1].real;
	Py[0] = abs2(tmp) / ((realnum) 2.0);

	Pz[0] = abs2(state[0]);
}

void spine::systems::system_1_singlet_triplet::updateHamiltonian(realnum * H)
{
	H[0] = exchange_interaction / ((realnum) 2.0);
	H[1] = magnetic_gradient / ((realnum) 2.0);
	H[2] = H[1];
	H[3] = -H[0];
}

#ifdef PLOT

void spine::systems::system_1_singlet_triplet::plotSetup(unsigned int points, double xmin, double xmax)
{
	this->points = points;
	this->xmin = xmin;
	this->xmax = xmax;
}

void spine::systems::system_1_singlet_triplet::plotAddU(complex * U, realnum t, bool plot_lab)
{
	unsigned int dim = getDimension();
	complex * statetmp = allocComplex(dim);
	complex * state = allocComplex(dim);

	// Update the overall operation
	initialize(state);
	initCopy(dim, statetmp, state);
	multiplyVector(dim, U, statetmp, state);
	plotAdd(state, t, plot_lab);

	freeComplex(state);
	freeComplex(statetmp);
}

void spine::systems::system_1_singlet_triplet::plotAdd(complex * state, realnum t, bool plot_lab)
{
	realnum Px, Py, Pz;
	measure(state, &Px, &Py, &Pz);

	if (p == nullptr)
	{
		p = new Plot(L"Qubit 0", points);
		p->add(255, 0, 0, xmin, xmax, 0, 1);
		p->add(0, 255, 0, xmin, xmax, 0, 1);
		p->add(0, 0, 255, xmin, xmax, 0, 1);
	}
	p->get(0)->add(t, Px);
	p->get(1)->add(t, Py);
	p->get(2)->add(t, Pz);
}

void spine::systems::system_1_singlet_triplet::plot()
{
	if (p != nullptr)
		p->redraw();
}

#endif

void spine::systems::system_1_singlet_triplet::setMagneticGradient(realnum value)
{
	magnetic_gradient = value;
}

realnum spine::systems::system_1_singlet_triplet::getMagneticGradient()
{
	return magnetic_gradient;
}

void spine::systems::system_1_singlet_triplet::setExchangeInteraction(realnum value)
{
	exchange_interaction = value;
}
