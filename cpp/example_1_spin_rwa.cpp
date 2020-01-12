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
#include "spine/solvers/solver_taylor.h"
#include "spine/systems/system_1_spin_rwa.h"

// Namespaces
using namespace std;
using namespace spine::math;

// Function prototypes
bool inHamiltonian(complex * H, realnum * dt);
void outOperation(complex * U);

// Private global variables/constants
static spine::systems::system_1_spin_rwa spin_system = spine::systems::system_1_spin_rwa();
static const realnum dt = (realnum) 50e-12;
static int Nsim, nsim;

int main(void)
{
	// Create the system
	spin_system.setMicrowaveFrequency((realnum)(2 * M_PI * 1e9 * 0));	// Set to zero for nicer plot
	spin_system.setLarmorFrequency((realnum)(2 * M_PI * 1e9 * 0));		// Set to zero for nicer plot
	spin_system.setRabiFrequency((realnum)(2 * M_PI * 1e6));

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
	spine::simulate(spin_system.getDimension(), inHamiltonian, outOperation, spine::solvers::solver_taylor);
	cout << "finished!" << endl;
#ifdef PLOT
	Plot::run();
#endif
}

bool inHamiltonian(complex * H, realnum * dt)
{
	if (nsim < Nsim)
	{
		// Set the signal at this time instance
		spin_system.setMicrowaveAmplitude(1);
		if (nsim < (Nsim / 2))
			spin_system.setMicrowavePhase(0);
		else
			spin_system.setMicrowavePhase((realnum)(M_PI / 2.0));

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
	spin_system.plotAddU(U, nsim * dt);
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
		// The ideal X, pi / 2 rotation
		complex Ux[4] = { {cos(M_PI / 4), 0}, {sin(0) * sin(M_PI / 4), -cos(0) * sin(M_PI / 4)}, {-sin(0) * sin(M_PI / 4), -cos(0) * sin(M_PI / 4)}, {cos(M_PI / 4), 0} };
		complex Uy[4] = { {cos(M_PI / 4), 0}, {sin(M_PI / 2) * sin(M_PI / 4), -cos(M_PI / 2) * sin(M_PI / 4)}, {-sin(M_PI / 2) * sin(M_PI / 4), -cos(M_PI / 2) * sin(M_PI / 4)}, {cos(M_PI / 4), 0} };
		complex Uideal[4];
		multiply(2, Uy, Ux, Uideal);

		// Last point, print the fidelity
		realnum F = spine::fidelity(spin_system.getDimension(), U, Uideal);
		cout << "Fidelity: " << F << endl;
	}
}
