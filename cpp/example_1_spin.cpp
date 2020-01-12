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
#include "spine/solvers/solver_analytical_xz.h"
#include "spine/systems/system_1_spin.h"

// Namespaces
using namespace std;
using namespace spine::math;

// Function prototypes
bool inHamiltonian(realnum * H, realnum * dt);
void outOperation(complex * U);

// Private global variables/constants
static spine::systems::system_1_spin spin_system = spine::systems::system_1_spin();
static const realnum dt = (realnum) 10e-12;
static int Nsim, nsim;

int main(void)
{
	// Create the system
	spin_system.setLarmorFrequency((realnum)(2 * M_PI * 1e9));
	spin_system.setRabiFrequency((realnum)(2 * M_PI * 10e6));

	// Generate the driving signal
	realnum tpi = (realnum)M_PI / spin_system.getRabiFrequency();
	int Npi = (int)(tpi / dt);

	// Simulate the system
	nsim = 0;
	Nsim = Npi;
#ifdef PLOT
	spin_system.plotSetup(Nsim, 0, tpi);
#endif
	cout << "Simulating..." << endl;
	spine::simulate(spin_system.getDimension(), inHamiltonian, outOperation, spine::solvers::solver_analytical_xz);
	cout << "finished!" << endl;
#ifdef PLOT
	Plot::run();
#endif
}

bool inHamiltonian(realnum * H, realnum * dt)
{
	if (nsim < Nsim)
	{
		// Set the signal at this time instance
		spin_system.setMicrowaveControl(cos(spin_system.getLarmorFrequency() * (nsim * ::dt - 0.5 * ::dt)));

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
	const unsigned int Nplot = 5;

	// Print the measurement probability
#ifdef PLOT
	spin_system.plotAddU(U, nsim * dt, true);
#endif
	if ((nsim % Nplot) == 0)
	{
		realnum P = abs2(U[0]);
		cout << P << endl;
#ifdef PLOT
		spin_system.plot();
#endif
	}

	// Continue the simulation
	nsim++;
	if (nsim == Nsim)
	{
		// Last point, print the fidelity
		realnum F = spine::fidelity(spin_system.getDimension(), U, (realnum) M_PI, (realnum) 0.0);
		cout << "Fidelity: " << F << endl;
	}
}
