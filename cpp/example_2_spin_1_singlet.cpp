/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

// Performance/accuracy settings
//#define MKL					(define at project level)
//#define DOUBLE_PRECISION		(define at project level)

// Includes
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include "spine/simulate.h"
#include "spine/fidelity.h"
#include "spine/math.h"
#include "spine/solvers/solver_diagonalization.h"
#include "spine/systems/system_2_spin_1_singlet.h"

// Namespaces
using namespace std;
using namespace spine::math;

// Function prototypes
bool inHamiltonian(realnum * H, realnum * dt);
void outOperation(complex * U);

// Private global variables/constants
static spine::systems::system_2_spin_1_singlet spin_system = spine::systems::system_2_spin_1_singlet();
static const realnum dt = (realnum) 50e-12;
static int Nsim, nsim;
static int Ngauss;
static double * ygauss;

int main(void)
{
	// Create the system
	spin_system.setChargingEnergy((realnum)(2 * M_PI * 100e9));
	spin_system.setRabiFrequency((realnum)(2 * M_PI * 1e6));
	spin_system.setTunnelControl((realnum)(2 * M_PI * 10e6 * sqrt(2)));
	spin_system.setLarmorFrequency(0, (realnum)(2 * M_PI * (1e9 - 10e6)));
	spin_system.setLarmorFrequency(1, (realnum)(2 * M_PI * (1e9 + 10e6)));

	// Generate the driving signal
	Ngauss = (int)(M_PI / 2 / spin_system.getRabiFrequency() / dt);
	ygauss = new double[Ngauss];
	double sum = 0;
	for (int n = 0; n < Ngauss; n++) {
		double x = -2.5 + 5.0 * n / (realnum)(Ngauss - 1);
		ygauss[n] = exp(-x * x / 2.0);
		sum += ygauss[n];
	}
	double scale = Ngauss / sum;
	for (int n = 0; n < Ngauss; n++)
		ygauss[n] = ygauss[n] * scale;

	// Simulate the system
	nsim = 0;
	Nsim = 7*Ngauss;
#ifdef PLOT
	spin_system.plotSetup(Nsim, 0, Nsim * dt);
#endif
	cout << "Simulating..." << endl;
	complex * state = allocComplex(spin_system.getDimension());
	spin_system.initialize(state);
	spine::simulate(spin_system.getDimension(), inHamiltonian, outOperation, spine::solvers::solver_diagonalization, state);
	freeComplex(state);
	cout << "finished!" << endl;
#ifdef PLOT
	Plot::run();
#endif

	delete[] ygauss;
}

bool inHamiltonian(realnum * H, realnum * dt)
{
	if (nsim < Nsim)
	{
		// Set the signal at this time instance
		spin_system.setMicrowaveControl(0);
		spin_system.setDetuningControl(0);
		if (nsim < Ngauss) {
			spin_system.setMicrowaveControl(0, ygauss[nsim] * cos(spin_system.getLarmorFrequency(0) * nsim * ::dt));
		}
		else if (nsim < 2 * Ngauss) {
			spin_system.setDetuningControl((1 - 2e-3) * spin_system.getChargingEnergy());
		}
		else if (nsim < 3 * Ngauss) {
			spin_system.setMicrowaveControl(0, ygauss[nsim - 2 * Ngauss] * cos(spin_system.getLarmorFrequency(0) * nsim * ::dt + M_PI / 2.0));
		}
		else if (nsim < 4 * Ngauss) {
			spin_system.setMicrowaveControl(1, 2 * ygauss[nsim - 3 * Ngauss] * cos(spin_system.getLarmorFrequency(1) * nsim * ::dt));
		}
		else if (nsim < 5 * Ngauss) {
			spin_system.setMicrowaveControl(0, ygauss[nsim - 4 * Ngauss] * cos(spin_system.getLarmorFrequency(0) * nsim * ::dt));
		}
		else if (nsim < 6 * Ngauss) {
			spin_system.setDetuningControl((1 - 2e-3) * spin_system.getChargingEnergy());
		}
		else if (nsim < 7 * Ngauss) {
			spin_system.setMicrowaveControl(0, ygauss[nsim - 6 * Ngauss] * cos(spin_system.getLarmorFrequency(0) * nsim * ::dt + M_PI / 2.0));
		}

		// Update the Hamiltonian accordingly
		spin_system.updateHamiltonian(H);

		// Provide the current timestep, and continue the simulation
		*dt = ::dt;
		return true;
	}
	return false;
}

void outOperation(complex * U)
{
	// Plot every Nplot-th point
	const unsigned int Nplot = 10;

	// Print the measurement probability
#ifdef PLOT
	spin_system.plotAdd(U, nsim * dt);
#endif
	if ((nsim % Nplot) == 0)
	{
		realnum P = abs2(U[4]);
		cout << P << endl;
#ifdef PLOT
		spin_system.plot();
#endif
	}

	// Continue the simulation
	nsim++;
}
