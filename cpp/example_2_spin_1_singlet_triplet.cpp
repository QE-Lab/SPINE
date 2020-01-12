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
#include "spine/systems/system_2_spin_1_singlet_triplet.h"

// Namespaces
using namespace std;
using namespace spine::math;

// Function prototypes
bool inHamiltonian(realnum * H, realnum * dt);
void outOperation(complex * U);

// Private global variables/constants
static spine::systems::system_2_spin_1_singlet_triplet spin_system = spine::systems::system_2_spin_1_singlet_triplet();
static const realnum dt = (realnum) 50e-12;
static int Nsim, nsim;
static int Ngauss;
static double * ygauss;
static double * ydetun;
static double * probability;

int main(void)
{
	// Create the system
	spin_system.setChargingEnergy((realnum)(2 * M_PI * 100e9));
	spin_system.setSingletTripletEnergy((realnum)(2 * M_PI * 10e9));
	spin_system.setRabiFrequency((realnum)(2 * M_PI * 1e6));
	spin_system.setTunnelControl((realnum)(2 * M_PI * 100e6 * sqrt(2)));
	spin_system.setLarmorFrequency(0, (realnum)(2 * M_PI * (1e9 - 10e6)));
	spin_system.setLarmorFrequency(1, (realnum)(2 * M_PI * (1e9 + 10e6)));

	// Generate the driving signal
	Ngauss = (int)(M_PI / spin_system.getRabiFrequency() / dt);
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

	double scale2 = (spin_system.getChargingEnergy() + spin_system.getSingletTripletEnergy() / 2);
	int Nadiabatic = (int)(M_PI / spin_system.getRabiFrequency() / 4 / dt);
	double * yadiabatic = new double[Nadiabatic];
	for (int n = 0; n < Nadiabatic; n++)
		yadiabatic[n] = (1.0 - exp(-n / ((double) Nadiabatic / 5.0))) / (1.0 - exp(-5.0)) * scale2;
	ydetun = new double[Ngauss];	
	for (int n = 0; n < Nadiabatic; n++)
		ydetun[n] = yadiabatic[n];
	for (int n = 0; n < Nadiabatic; n++)
		ydetun[Ngauss - 1 - n] = yadiabatic[n];
	for (int n = 0; n < Ngauss - 2 * Nadiabatic; n++)
		ydetun[Nadiabatic + n] = scale2;
	delete[] yadiabatic;

	// Simulate the system
	nsim = 0;
	Nsim = 3*Ngauss;
	probability = new double[Nsim];
#ifdef PLOT
	Plot(L"Gaussian", Ngauss, ygauss, 0, 0, 0, 0, Ngauss, 0, scale);
	Plot(L"Detuning", Ngauss, ydetun, 0, 0, 0, 0, Ngauss, 0, scale2);
	spin_system.plotSetup(Nsim, 0, Nsim * dt);
#endif
	cout << "Simulating..." << endl;
	spine::simulate(spin_system.getDimension(), inHamiltonian, outOperation, spine::solvers::solver_diagonalization);
	cout << "finished!" << endl;
#ifdef PLOT
	Plot(L"Probability", Nsim, probability, 0, 0, 0, 0, Nsim, 0, 1);
	Plot::run();
#endif

	delete[] ygauss;
	delete[] ydetun;
	delete[] probability;
}

bool inHamiltonian(realnum * H, realnum * dt)
{
	if (nsim < Nsim)
	{
		// Set the signal at this time instance
		spin_system.setMicrowaveControl(0);
		spin_system.setDetuningControl(0);
		if (nsim < Ngauss) {
			spin_system.setDetuningControl(ydetun[nsim]);
		}
		else if (nsim < 2 * Ngauss) {
			spin_system.setMicrowaveControl(0, ygauss[nsim - 1 * Ngauss] * cos(spin_system.getLarmorFrequency(0) * nsim * ::dt));
		}
		else if (nsim < 3 * Ngauss) {
			spin_system.setDetuningControl(ydetun[nsim - 2 * Ngauss]);
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
	probability[nsim] = abs2(U[8*4]) + abs2(U[8*5]) + abs2(U[8*6]) + abs2(U[8*7]);
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
