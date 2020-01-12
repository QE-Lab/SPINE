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
#include "spine/systems/system_n_spin_n_singlet.h"

// Namespaces
using namespace std;
using namespace spine::math;

// Function prototypes
bool inHamiltonian(realnum * H, realnum * dt);
void outOperation(complex * U);

// Private global variables/constants
static spine::systems::system_n_spin_n_singlet spin_system = spine::systems::system_n_spin_n_singlet(4);
static const realnum dt = (realnum) 50e-12;
static int Nsim, nsim;
static int Ngauss;
static realnum * microwave0;
static realnum * microwave1;
static realnum * microwave2;
static realnum * microwave3;
static realnum * tunnel01;
static realnum * tunnel02;
static realnum * tunnel13;
static realnum * tunnel23;
static realnum * detuning2;
static realnum * detuning3;

int main(void)
{
	// Set the charging energy of the dots
	spin_system.setChargingEnergy(0, 2 * M_PI * 100e9);
	spin_system.setChargingEnergy(1, 2 * M_PI * 100e9);
	spin_system.setChargingEnergy(2, 2 * M_PI * 100e9);
	spin_system.setChargingEnergy(3, 2 * M_PI * 100e9);

	// Set the Larmor frequencies of the dots
	spin_system.setLarmorFrequency(0, 2 * M_PI * 2.5e9);
	spin_system.setLarmorFrequency(1, 2 * M_PI * 2.5e9);
	spin_system.setLarmorFrequency(2, 2 * M_PI * (2.5e9 + 0.5e9*sqrt(2.0)));
	spin_system.setLarmorFrequency(3, 2 * M_PI * (2.5e9 - 0.5e9*sqrt(2.0)));

	// Tunnel rates for different operations
	realnum tunnel_coupling_swap = 2 * M_PI * 0.5e9;
	realnum tunnel_coupling_cz = 2 * M_PI * 0.5e9 * 1.5;
	realnum tunnel_coupling_measure = 2 * M_PI * 0.5e9 * 2;

	// The system is tuned such that each operation has a speed of 10 MHz
	realnum omega_SWAP = 2 * M_PI * 10e6;
	realnum omega_CZ = 2 * M_PI * 10e6;
	realnum omega_R = 2 * M_PI * 10e6;
	spin_system.setRabiFrequency(omega_R);

	// Perform the following simulation
	unsigned int Npi = (unsigned int)round(M_PI / omega_R / dt);
	unsigned int N0 = 0;             // -			 initialize the system to the ground state
	unsigned int N1 = N0 + Npi;      // 50.0 ns      rotate q1 from 0 to 1 using a microwave signal
	unsigned int N2 = N1 + Npi;      // 50.0 ns      swap q1 and q2 using the tunnel coupling
	unsigned int N3 = N2 + Npi / 2;  // 25.0 ns      rotate q3 and q4 in plane using a microwave signal
	unsigned int N4 = N3 + Npi / 4;  // 12.5 ns      adiabatic change the tunnel coupling between q1 / q3 and q2 / q4
	unsigned int N5 = N4 + Npi / 2;  // 25.0 ns      controlled - z between q1 / q3 and q2 / q4
	unsigned int N6 = N5 + Npi / 4;  // 12.5 ns      adiabatic change the tunnel coupling between q1 / q3 and q2 / q4
	unsigned int N7 = N6 + Npi / 2;  // 25.0 ns      rotate q3 and q4 out plane using a microwave signal, out of phase
	unsigned int N8 = N7 + Npi;      // 50.0 ns      measure q3 / q4
	unsigned int N9 = N8 + Npi;      // 50.0 ns      rotate q1 using a microwave signal
	unsigned int N = N9;

	// Generate the I / Q clocks for the simulation
	realnum * t = new realnum[N];
	realnum * clk_I0 = new realnum[N];
	realnum * clk_Q0 = new realnum[N];
	realnum * clk_I1 = new realnum[N];
	realnum * clk_Q1 = new realnum[N];
	realnum * clk_I2 = new realnum[N];
	realnum * clk_Q2 = new realnum[N];
	realnum * clk_I3 = new realnum[N];
	realnum * clk_Q3 = new realnum[N];
	for (unsigned int i = 0; i < N; i++) {
		t[i] = i * dt;
		clk_I0[i] = sin(spin_system.getLarmorFrequency(0) * t[i]);
		clk_Q0[i] = cos(spin_system.getLarmorFrequency(0) * t[i]);
		clk_I1[i] = sin(spin_system.getLarmorFrequency(1) * t[i]);
		clk_Q1[i] = cos(spin_system.getLarmorFrequency(1) * t[i]);
		clk_I2[i] = sin(spin_system.getLarmorFrequency(2) * t[i]);
		clk_Q2[i] = cos(spin_system.getLarmorFrequency(2) * t[i]);
		clk_I3[i] = sin(spin_system.getLarmorFrequency(3) * t[i]);
		clk_Q3[i] = cos(spin_system.getLarmorFrequency(3) * t[i]);
	}

	// All signals
	microwave0 = new realnum[N]();
	microwave1 = new realnum[N]();
	microwave2 = new realnum[N]();
	microwave3 = new realnum[N]();
	tunnel01 = new realnum[N]();
	tunnel02 = new realnum[N]();
	tunnel13 = new realnum[N]();
	tunnel23 = new realnum[N]();
	detuning2 = new realnum[N]();
	detuning3 = new realnum[N]();

	// Apply the microwave signal to q0
	for (unsigned int n = N0; n < N1; n++)
		microwave0[n] = clk_I0[n];

	// Change the tunnel coupling diabatically for a SWAP gate
	for (unsigned int n = N1; n < N2; n++)
		tunnel01[n] = tunnel_coupling_swap;

	// Apply the microwave signal to q2 and q3
	for (unsigned int n = N2; n < N3; n++) {
		microwave2[n] = clk_I2[n];
		microwave3[n] = clk_I3[n];
	}

	// Ramp the tunnel coupling up for the CZ
	for (unsigned int n = N3; n < N4; n++) {
		realnum t0 = ((realnum)n - (realnum)N3) / ((realnum)N4 - (realnum)N3) * tunnel_coupling_cz;
		tunnel02[n] = t0;
		tunnel13[n] = t0;
		detuning2[n] = 2 * M_PI * 75e9;
		detuning3[n] = 2 * M_PI * 75e9;
	}

	//  Hold the tunnel coupling for the CZ
	for (unsigned int n = N4; n < N5; n++) {
		tunnel02[n] = tunnel_coupling_cz;
		tunnel13[n] = tunnel_coupling_cz;
		detuning2[n] = 2 * M_PI * 75e9;
		detuning3[n] = 2 * M_PI * 75e9;
	}

	// Ramp the tunnel coupling down for the CZ
	for (unsigned int n = N5; n < N6; n++) {
		realnum t0 = (1.0 - ((realnum)n - (realnum)N5) / ((realnum)N6 - (realnum)N5)) * tunnel_coupling_cz;
		tunnel02[n] = t0;
		tunnel13[n] = t0;
		detuning2[n] = 2 * M_PI * 75e9;
		detuning3[n] = 2 * M_PI * 75e9;
	}

	// Apply the microwave signal to q2 and q3
	for (unsigned int n = N6; n < N7; n++) {
		microwave2[n] = clk_Q2[n];
		microwave3[n] = clk_Q3[n];
		detuning2[n] = 2 * M_PI * 75e9;
		detuning3[n] = 2 * M_PI * 75e9;
	}

	// Detune and ramp the tunnel coupling
	for (unsigned int n = N7; n < N8; n++) {
		realnum rel = ((realnum)n - (realnum)N7) / ((realnum)N8 - (realnum)N7);
		detuning2[n] = 2 * M_PI * 75e9 * (1.0 - rel);
		detuning3[n] = 2 * M_PI * 75e9 * (1.0 + rel);
		tunnel23[n] = rel * tunnel_coupling_measure;
	}

	// Apply the microwave signal to q1
	for (unsigned int n = N8; n < N9; n++)
		microwave1[n] = clk_Q1[n];

	// Simulate the system
	nsim = 0;
	Nsim = N;
#ifdef PLOT
	spin_system.plotSetup(N, 0, N * dt);
#endif
	cout << "Simulating..." << endl;
	spine::simulate(spin_system.getDimension(), inHamiltonian, outOperation, spine::solvers::solver_diagonalization);
	cout << "finished!" << endl;
#ifdef PLOT
	Plot::run();
#endif

	// Free the allocated memory
	delete[] t;
	delete[] clk_I0;
	delete[] clk_Q0;
	delete[] clk_I1;
	delete[] clk_Q1;
	delete[] clk_I2;
	delete[] clk_Q2;
	delete[] clk_I3;
	delete[] clk_Q3;
	delete[] microwave0;
	delete[] microwave1;
	delete[] microwave2;
	delete[] microwave3;
	delete[] tunnel01;
	delete[] tunnel02;
	delete[] tunnel13;
	delete[] tunnel23;
	delete[] detuning2;
	delete[] detuning3;
}

bool inHamiltonian(realnum * H, realnum * dt)
{
	if (nsim < Nsim)
	{
		// Set the signal at this time instance
		spin_system.setMicrowaveControl(0, microwave0[nsim]);
		spin_system.setMicrowaveControl(1, microwave1[nsim]);
		spin_system.setMicrowaveControl(2, microwave2[nsim]);
		spin_system.setMicrowaveControl(3, microwave3[nsim]);
		spin_system.setTunnelControl(0, 1, tunnel01[nsim]);
		spin_system.setTunnelControl(0, 2, tunnel02[nsim]);
		spin_system.setTunnelControl(1, 3, tunnel13[nsim]);
		spin_system.setTunnelControl(2, 3, tunnel23[nsim]);
		spin_system.setDetuningControl(2, detuning2[nsim]);
		spin_system.setDetuningControl(3, detuning3[nsim]);
		
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
	spin_system.plotAddU(U, nsim * dt);
	if ((nsim % Nplot) == 0)
		spin_system.plot();
#endif

	// Continue the simulation
	nsim++;
}
