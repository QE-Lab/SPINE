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
#include "spine/systems/system_1_singlet_triplet.h"

// Namespaces
using namespace std;
using namespace spine::math;

// Function prototypes
bool inHamiltonian(realnum * H, realnum * dt);
void outOperation(complex * U);

// Private global variables/constants
static spine::systems::system_1_singlet_triplet spin_system = spine::systems::system_1_singlet_triplet();
static int Nsim, nsim;
static realnum J[2];
static realnum t[2];

int main(void)
{
	// Create the system
	spin_system.setMagneticGradient(160e6);

	// Generate the driving signal, consisting of two levels with different duration
	J[0] = spin_system.getMagneticGradient() * (1.0 + sqrt(2.0));
	J[1] = spin_system.getMagneticGradient() * (J[0] - spin_system.getMagneticGradient()) / (J[0] + spin_system.getMagneticGradient());
	t[0] = M_PI / sqrt(spin_system.getMagneticGradient() * spin_system.getMagneticGradient() + J[0] * J[0]);
	t[1] = M_PI / sqrt(spin_system.getMagneticGradient() * spin_system.getMagneticGradient() + J[1] * J[1]);

	// Simulate the system
	nsim = 0;
	Nsim = 50;
#ifdef PLOT
	spin_system.plotSetup(Nsim, 0, t[0] + t[1]);
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
		spin_system.setExchangeInteraction(J[0]);

		// Update the Hamiltonian accordingly
		spin_system.updateHamiltonian(H);

		// Provide the current timestep, and continue the simulation
		*dt = t[0] / Nsim;
		return true;
	}
	else if (nsim < 2 * Nsim)
	{
		// Set the signal at this time instance
		spin_system.setExchangeInteraction(J[1]);

		// Update the Hamiltonian accordingly
		spin_system.updateHamiltonian(H);

		// Provide the current timestep, and continue the simulation
		*dt = t[1] / Nsim;
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
	realnum time;
	if (nsim < Nsim)
		time = t[0] / Nsim * nsim;
	else
		time = t[0] + t[1] / Nsim * (nsim - Nsim);
	spin_system.plotAddU(U, time);
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
	if (nsim == 2*Nsim)
	{
		// Last point, print the fidelity
		realnum F = spine::fidelity(spin_system.getDimension(), U, (realnum)(M_PI / 2.0), (realnum)(-M_PI / 2.0));
		cout << "Fidelity: " << F << endl;
	}
}
