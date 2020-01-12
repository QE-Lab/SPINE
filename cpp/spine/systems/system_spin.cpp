/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "system_spin.h"

using namespace spine;

spine::systems::system_spin::system_spin(unsigned int dots)
{
	this->dots = dots;
	this->larmor_frequency = new realnum[dots]();
	this->rabi_frequency = new realnum[dots]();
	this->charging_energy = new realnum[dots]();
	this->singlet_triplet_energy = new realnum[dots]();
	this->microwave_control = new realnum[dots]();
	this->detuning_control = new realnum[dots]();
	this->tunnel_control = new realnum[dots*dots]();
}

spine::systems::system_spin::~system_spin()
{
	delete[] larmor_frequency;
	delete[] rabi_frequency;
	delete[] charging_energy;
	delete[] singlet_triplet_energy;
	delete[] microwave_control;
	delete[] detuning_control;
	delete[] tunnel_control;
#ifdef PLOT
	if (p != nullptr) {
		for (unsigned int dot = 0; dot < dots; dot++)
			delete p[dot];
		delete[] p;
	}
#endif
}

int spine::systems::system_spin::initializeIndex()
{
	unsigned int * states = new unsigned int[dots];
	for (unsigned int dot = 0; dot < dots; dot++)
		states[dot] = STATE_0;
	int index = getIndex(states);
	delete[] states;
	return index;
}

void spine::systems::system_spin::initialize(complex * state)
{
	// Zero the state vector
	initZero(getDimension(), state);

	// Set the state at the index of the ground state to 1
	int index = initializeIndex();
	if (index >= 0)
		state[index].real = (realnum) 1.0;
}

void spine::systems::system_spin::measure(complex * state, realnum t, realnum * Px, realnum * Py, realnum * Pz)
{
	// For every dot
	for (unsigned int dot = 0; dot < dots; dot++) {
		Px[dot] = 0;
		Py[dot] = 0;
		Pz[dot] = 0;

		// Measure in rotating frame at time t if requested
		complex rot_x = { 1, 0 };
		complex rot_y = { 0, 1 };
		if (t != 0) {
#ifdef DOUBLE_PRECISION
			rot_x.real = cos(larmor_frequency[dot] * t);
			rot_x.imag = sin(larmor_frequency[dot] * t);
			rot_y.real = cos(larmor_frequency[dot] * t + M_PI / 2.0);
			rot_y.imag = sin(larmor_frequency[dot] * t + M_PI / 2.0);
#else
			rot_x.real = cosf(larmor_frequency[dot] * t);
			rot_x.imag = sinf(larmor_frequency[dot] * t);
			rot_y.real = cosf(larmor_frequency[dot] * t + (realnum)M_PI / 2.0f);
			rot_y.imag = sinf(larmor_frequency[dot] * t + (realnum)M_PI / 2.0f);
#endif
		}

		// Sum over all indices where this dot in the requested state
		vector<unsigned int> * indices0 = getIndexMeasurement(dot, STATE_0);
		vector<unsigned int> * indices1 = getIndexMeasurement(dot, STATE_1);
		vector<unsigned int>::iterator it1 = indices1->begin();
		for (vector<unsigned int>::iterator it0 = indices0->begin(); it0 != indices0->end(); ++it0) {
			unsigned int index0 = *it0;
			unsigned int index1 = *it1;
			complex tmp;

			// |(state[index0] + rot_x * state[index1]) / sqrt(2)| ^ 2
			tmp.real = state[index0].real + rot_x.real * state[index1].real - rot_x.imag * state[index1].imag;
			tmp.imag = state[index0].imag + rot_x.imag * state[index1].real + rot_x.real * state[index1].imag;
			Px[dot] += abs2(tmp) / ((realnum) 2.0);

			// |(state[index0] + rot_y * state[index1]) / sqrt(2)| ^ 2
			tmp.real = state[index0].real + rot_y.real * state[index1].real - rot_y.imag * state[index1].imag;
			tmp.imag = state[index0].imag + rot_y.imag * state[index1].real + rot_y.real * state[index1].imag;
			Py[dot] += abs2(tmp) / ((realnum) 2.0);

			// |state[index0]| ^ 2
			Pz[dot] += abs2(state[index0]);

			++it1;
		}
	}
}

void spine::systems::system_spin::measure(complex * state, realnum * Px, realnum * Py, realnum * Pz)
{
	return measure(state, 0, Px, Py, Pz);
}

void spine::systems::system_spin::measureST(complex * state, realnum * Pn, realnum * Ps, realnum * Pt)
{
	// For every dot
	for (unsigned int dot = 0; dot < dots; dot++) {
		Pn[dot] = 0;
		Ps[dot] = 0;
		Pt[dot] = 0;

		// Sum over all indices where this dot in the requested state
		vector<unsigned int> * indices = getIndexMeasurement(dot, STATE_N);
		for (vector<unsigned int>::iterator it = indices->begin(); it != indices->end(); ++it)
			Pn[dot] += abs2(state[*it]);
		indices = getIndexMeasurement(dot, STATE_S);
		for (vector<unsigned int>::iterator it = indices->begin(); it != indices->end(); ++it)
			Ps[dot] += abs2(state[*it]);
		indices = getIndexMeasurement(dot, STATE_T0);
		for (vector<unsigned int>::iterator it = indices->begin(); it != indices->end(); ++it)
			Pt[dot] += abs2(state[*it]);
		indices = getIndexMeasurement(dot, STATE_TP);
		for (vector<unsigned int>::iterator it = indices->begin(); it != indices->end(); ++it)
			Pt[dot] += abs2(state[*it]);
		indices = getIndexMeasurement(dot, STATE_TM);
		for (vector<unsigned int>::iterator it = indices->begin(); it != indices->end(); ++it)
			Pt[dot] += abs2(state[*it]);
	}
}

#ifdef PLOT

void spine::systems::system_spin::plotSetup(unsigned int points, double xmin, double xmax)
{
	this->points = points;
	this->xmin = xmin;
	this->xmax = xmax;
}

void spine::systems::system_spin::plotAddU(complex * U, realnum t, bool plot_lab)
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

void spine::systems::system_spin::plotAdd(complex * state, realnum t, bool plot_lab)
{
	realnum * Px = nullptr;
	realnum * Py = nullptr;
	realnum * Pz = nullptr;
	if (plot_lab)
	{
		Px = new realnum[dots];
		Py = new realnum[dots];
		Pz = new realnum[dots];
		measure(state, Px, Py, Pz);
	}
	realnum * Pxr = new realnum[dots];
	realnum * Pyr = new realnum[dots];
	realnum * Pzr = new realnum[dots];
	measure(state, t, Pxr, Pyr, Pzr);

	if (p == nullptr)
	{
		p = new Plot*[dots];
		for (unsigned int dot = 0; dot < dots; dot++)
		{
			wchar_t title[10];
			swprintf(title, 10, L"Qubit %d", dot);
			p[dot] = new Plot(title, points);
			if (plot_lab)
			{
				p[dot]->add(255, 192, 192, xmin, xmax, 0, 1);
				p[dot]->add(192, 255, 192, xmin, xmax, 0, 1);
			}
			p[dot]->add(255, 0, 0, xmin, xmax, 0, 1);
			p[dot]->add(0, 255, 0, xmin, xmax, 0, 1);
			p[dot]->add(0, 0, 255, xmin, xmax, 0, 1);
		}
	}
	for (unsigned int dot = 0; dot < dots; dot++)
	{
		if (plot_lab)
		{
			p[dot]->get(0)->add(t, Px[dot]);
			p[dot]->get(1)->add(t, Py[dot]);
			p[dot]->get(2)->add(t, Pxr[dot]);
			p[dot]->get(3)->add(t, Pyr[dot]);
			p[dot]->get(4)->add(t, Pzr[dot]);
		}
		else
		{
			p[dot]->get(0)->add(t, Pxr[dot]);
			p[dot]->get(1)->add(t, Pyr[dot]);
			p[dot]->get(2)->add(t, Pzr[dot]);
		}
	}

	if (plot_lab)
	{
		delete[] Px;
		delete[] Py;
		delete[] Pz;
	}
	delete[] Pxr;
	delete[] Pyr;
	delete[] Pzr;
}

void spine::systems::system_spin::plot()
{
	if (p != nullptr)
		for (unsigned int dot = 0; dot < dots; dot++)
			p[dot]->redraw();
}

#endif

void spine::systems::system_spin::setLarmorFrequency(realnum value)
{
	for (unsigned int dot = 0; dot < dots; dot++)
		larmor_frequency[dot] = value;
}

void spine::systems::system_spin::setLarmorFrequency(unsigned int dot, realnum value)
{
	larmor_frequency[dot] = value;
}

realnum spine::systems::system_spin::getLarmorFrequency()
{
	return larmor_frequency[0];
}

realnum spine::systems::system_spin::getLarmorFrequency(unsigned int dot)
{
	return larmor_frequency[dot];
}

void spine::systems::system_spin::setRabiFrequency(realnum value)
{
	for (unsigned int dot = 0; dot < dots; dot++)
		rabi_frequency[dot] = value;
}

void spine::systems::system_spin::setRabiFrequency(unsigned int dot, realnum value)
{
	rabi_frequency[dot] = value;
}

realnum spine::systems::system_spin::getRabiFrequency()
{
	return rabi_frequency[0];
}

realnum spine::systems::system_spin::getRabiFrequency(unsigned int dot)
{
	return rabi_frequency[dot];
}

void spine::systems::system_spin::setChargingEnergy(realnum value)
{
	for (unsigned int dot = 0; dot < dots; dot++)
		charging_energy[dot] = value;
}

void spine::systems::system_spin::setChargingEnergy(unsigned int dot, realnum value)
{
	charging_energy[dot] = value;
}

realnum spine::systems::system_spin::getChargingEnergy()
{
	return charging_energy[0];
}

realnum spine::systems::system_spin::getChargingEnergy(unsigned int dot)
{
	return charging_energy[dot];
}

void spine::systems::system_spin::setSingletTripletEnergy(realnum value)
{
	for (unsigned int dot = 0; dot < dots; dot++)
		singlet_triplet_energy[dot] = value;
}

void spine::systems::system_spin::setSingletTripletEnergy(unsigned int dot, realnum value)
{
	singlet_triplet_energy[dot] = value;
}

realnum spine::systems::system_spin::getSingletTripletEnergy()
{
	return singlet_triplet_energy[0];
}

realnum spine::systems::system_spin::getSingletTripletEnergy(unsigned int dot)
{
	return singlet_triplet_energy[dot];
}

void spine::systems::system_spin::setMicrowaveControl(realnum value)
{
	for (unsigned int dot = 0; dot < dots; dot++)
		microwave_control[dot] = value;
}

void spine::systems::system_spin::setMicrowaveControl(unsigned int dot, realnum value)
{
	microwave_control[dot] = value;
}

void spine::systems::system_spin::setDetuningControl(realnum value)
{
	for (unsigned int dot = 0; dot < dots; dot++)
		detuning_control[dot] = value;
}

void spine::systems::system_spin::setDetuningControl(unsigned int dot, realnum value)
{
	detuning_control[dot] = value;
}

void spine::systems::system_spin::setTunnelControl(realnum value)
{
	for (unsigned int dota = 0; dota < dots; dota++) {
		for (unsigned int dotb = 0; dotb < dots; dotb++) {
			tunnel_control[dota + dots * dotb] = value;
			tunnel_control[dotb + dots * dota] = value;
		}
	}
}

void spine::systems::system_spin::setTunnelControl(unsigned int dota, unsigned int dotb, realnum value)
{
	tunnel_control[dota + dots * dotb] = value;
	tunnel_control[dotb + dots * dota] = value;
}
