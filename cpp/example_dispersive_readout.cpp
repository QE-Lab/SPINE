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
#include "spine/systems/system_dispersive_readout.h"

// Namespaces
using namespace std;
using namespace spine::math;

// Function prototypes
bool inHamiltonian(complex * H, realnum * dt);
void outOperation(complex * U);

// Private global variables/constants
static spine::systems::system_dispersive_readout spin_system = spine::systems::system_dispersive_readout();
static int Nsim, nsim;
static realnum dt;
static realnum * gate_voltage;
static realnum * probability;

int main(void)
{
	// Create the system
	spin_system.setElectronTemperature(0.1);
	spin_system.setTunnelRate((realnum)(2.0 * M_PI * 1e9));
	spin_system.setLeverArm(0.05);

	// Generate the driving signal
	// Drive frequency determines waveform shape(capacitive vs resistive regimes)
	realnum w0 = 2.0 * spin_system.getTunnelRate();
	realnum tmax = 10.0 * M_PI / w0;
	unsigned int Npts = 1000;
	dt = tmax / Npts;
	gate_voltage = new realnum[Npts];
	probability = new realnum[Npts];
	for (unsigned int i = 0; i < Npts; i++)
		gate_voltage[i]	= 0.4e-3 / spin_system.getLeverArm() * sin(w0 * dt * i);

	// Simulate the system
	nsim = 0;
	Nsim = Npts;
#ifdef PLOT
	spin_system.plotSetup(Nsim, 0, Nsim * dt);
#endif
	cout << "Simulating..." << endl;
	spine::simulate(spin_system.getDimension(), inHamiltonian, outOperation, spine::solvers::solver_taylor);
	cout << "finished!" << endl;

	// Plot the voltage and current
#ifdef PLOT
	realnum q = 1.602e-19;
	realnum * current = new realnum[Npts - 1];
	for (unsigned int i = 0; i < Npts - 1; i++)
		current[i] = -q * (probability[i + 1] - probability[i]) / dt;
	Plot(L"Voltage", Npts, gate_voltage, 0, 255, 0, 0, Nsim * dt, -0.01, 0.01);
	Plot(L"Current", Npts - 1, current, 255, 0, 0, 0, Nsim * dt, -1e-9, 1e-9);
	Plot::run();
#endif

	delete[] gate_voltage;
	delete[] probability;
	delete[] current;
}

bool inHamiltonian(complex * H, realnum * dt)
{
	if (nsim < Nsim)
	{
		// Set the signal at this time instance
		spin_system.setGateVoltage(gate_voltage[nsim]);

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
	probability[nsim] = U[2].real;

#ifdef PLOT
	spin_system.plotAddU(U, nsim * dt);
#endif
	if ((nsim % Nplot) == 0)
	{
		cout << probability[nsim] << endl;
#ifdef PLOT
		spin_system.plot();
#endif
	}

	// Continue the simulation
	nsim++;
}
