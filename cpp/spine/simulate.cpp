/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "simulate.h"

using namespace spine;

void spine::simulate(unsigned int dim,
	bool(*inHamiltonian)(realnum * H, realnum * timestep),
	void(*outOperation)(complex * U),
	void(*solver)(unsigned int dim, realnum * H, complex * dU))
{
	unsigned int N = dim * dim;

	// Allocate space for the operation
	complex * dU = allocComplex(N);
	complex * Utmp = allocComplex(N);
	complex * U = allocComplex(N);
	initIdentity(dim, U);

	// Allocate space for the Hamiltonian and timestep
	realnum * dtH = allocReal(N);
	realnum timestep;

	// Let the user set the Hamiltonian, actually (dt * H), the argument of expm()
	while (inHamiltonian(dtH, &timestep))
	{
		// Scale with the timestep
		for (unsigned int i = 0; i < N; i++)
			dtH[i] *= timestep;

		// Solve for the incremental operation
		solver(dim, dtH, dU);

		// Update the overall operation
		initCopy(N, Utmp, U);
		multiply(dim, dU, Utmp, U);

		// Callback with the resulting operation
		outOperation(U);
	}

	// Free the allocated memory
	freeReal(dtH);
	freeComplex(U);
	freeComplex(Utmp);
	freeComplex(dU);
}

void spine::simulate(unsigned int dim,
	bool(*inHamiltonian)(realnum * H, realnum * timestep),
	void(*outOperation)(complex * U),
	void(*solver)(unsigned int dim, realnum * H, complex * dU),
	complex * state)
{
	unsigned int N = dim * dim;

	// Allocate space for the operation
	complex * dU = allocComplex(N);
	complex * statetmp = allocComplex(dim);

	// Allocate space for the Hamiltonian and timestep
	realnum * dtH = allocReal(N);
	realnum timestep;

	// Let the user set the Hamiltonian, actually (dt * H), the argument of expm()
	while (inHamiltonian(dtH, &timestep))
	{
		// Scale with the timestep
		for (unsigned int i = 0; i < N; i++)
			dtH[i] *= timestep;

		// Solve for the incremental operation
		solver(dim, dtH, dU);

		// Update the overall operation
		initCopy(dim, statetmp, state);
		multiplyVector(dim, dU, statetmp, state);

		// Callback with the resulting operation
		outOperation(state);
	}

	// Free the allocated memory
	freeReal(dtH);
	freeComplex(statetmp);
	freeComplex(dU);
}

void spine::simulate(unsigned int dim,
	bool(*inHamiltonian)(complex * H, realnum * timestep),
	void(*outOperation)(complex * U),
	void(*solver)(unsigned int dim, complex * H, complex * dU))
{
	unsigned int N = dim * dim;

	// Allocate space for the operation
	complex * dU = allocComplex(N);
	complex * Utmp = allocComplex(N);
	complex * U = allocComplex(N);
	initIdentity(dim, U);

	// Allocate space for the Hamiltonian and timestep
	complex * dtH = allocComplex(N);
	realnum timestep;

	// Let the user set the Hamiltonian, actually (dt * H), the argument of expm()
	while (inHamiltonian(dtH, &timestep))
	{
		// Scale with the timestep
		for (unsigned int i = 0; i < N; i++) {
			dtH[i].real *= timestep;
			dtH[i].imag *= timestep;
		}

		// Solve for the incremental operation
		solver(dim, dtH, dU);

		// Update the overall operation
		initCopy(N, Utmp, U);
		multiply(dim, dU, Utmp, U);

		// Callback with the resulting operation
		outOperation(U);
	}

	// Free the allocated memory
	freeComplex(dtH);
	freeComplex(U);
	freeComplex(Utmp);
	freeComplex(dU);
}

void spine::simulate(unsigned int dim,
	bool(*inHamiltonian)(complex * H, realnum * timestep),
	void(*outOperation)(complex * U),
	void(*solver)(unsigned int dim, complex * H, complex * dU),
	complex * state)
{
	unsigned int N = dim * dim;

	// Allocate space for the operation
	complex * dU = allocComplex(N);
	complex * statetmp = allocComplex(dim);

	// Allocate space for the Hamiltonian and timestep
	complex * dtH = allocComplex(N);
	realnum timestep;

	// Let the user set the Hamiltonian, actually (dt * H), the argument of expm()
	while (inHamiltonian(dtH, &timestep))
	{
		// Scale with the timestep
		for (unsigned int i = 0; i < N; i++) {
			dtH[i].real *= timestep;
			dtH[i].imag *= timestep;
		}

		// Solve for the incremental operation
		solver(dim, dtH, dU);

		// Update the overall operation
		initCopy(dim, statetmp, state);
		multiplyVector(dim, dU, statetmp, state);

		// Callback with the resulting operation
		outOperation(state);
	}

	// Free the allocated memory
	freeComplex(dtH);
	freeComplex(statetmp);
	freeComplex(dU);
}
