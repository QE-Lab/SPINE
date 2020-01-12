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
#include "spine/systems/system_2_spin_2_singlet.h"

// Namespaces
using namespace std;
using namespace spine::math;

// Function prototypes
bool inHamiltonian(realnum * H, realnum * dt);
void outOperation(complex * U);

// Private global variables/constants
static spine::systems::system_2_spin_2_singlet spin_system = spine::systems::system_2_spin_2_singlet();
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
	spin_system.setLarmorFrequency(0, (realnum)(2 * M_PI * 1e9));
	spin_system.setLarmorFrequency(1, (realnum)(2 * M_PI * 1e9));

	// Generate the driving signal
	Nsim = (int)(500e-9 / dt);
	spin_system.setDetuningControl(1, (1 - 2e-3) * spin_system.getChargingEnergy());

	// Simulate the system
	nsim = 0;
#ifdef PLOT
	spin_system.plotSetup(Nsim, 0, Nsim * dt);
#endif
	cout << "Simulating..." << endl;
	complex * state = allocComplex(spin_system.getDimension());
	unsigned int states[2] = { spin_system.STATE_0, spin_system.STATE_0 };
	unsigned int index00 = spin_system.getIndex(states);
	states[0] = spin_system.STATE_1;
	unsigned int index01 = spin_system.getIndex(states);
	state[index00].real = 1.0 / sqrt(2.0);
	state[index01].real = 1.0 / sqrt(2.0);
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
		realnum P = abs2(U[4]) + abs2(U[5]);
		cout << P << endl;
#ifdef PLOT
		spin_system.plot();
#endif
	}

	// Continue the simulation
	nsim++;
}
