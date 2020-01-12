/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "system_n_spin_n_singlet.h"

using namespace spine;

spine::systems::system_n_spin_n_singlet::system_n_spin_n_singlet(unsigned int dots)
	: spine::systems::system_spin(dots)
{
	// Determine the valid codes
	size = (1 << (dots << 1));
	storage_lookup = new unsigned int[size];
	storage_size = 0;
	for (unsigned int i = 0; i < size; i++) {
		if (isValid(i)) {
			storage_lookup[i] = storage_size;
			storage_size++;
		}
		else {
			storage_lookup[i] = -1;
		}
	}

	// Determine where the measurements should take place
	probability_i = new vector<unsigned int>[dots * 4];

	// Determine where the values are to be places in the Hamiltonian
	charging_energies_i = new vector<unsigned int>[dots];
	larmor_frequencies_ip = new vector<unsigned int>[dots];
	larmor_frequencies_im = new vector<unsigned int>[dots];
	tunnel_couplings_ip = new vector<storage_location>[dots*dots];
	tunnel_couplings_im = new vector<storage_location>[dots*dots];
	microwave_control_i = new vector<storage_location>[dots];
	detuning_control_ip = new vector<unsigned int>[dots];
	detuning_control_im = new vector<unsigned int>[dots];

	for (unsigned int dot = 0; dot < dots; dot++)
	{
		// Create the bitmasks
		unsigned int is00 = 0x00 * (1 << (dot << 1));
		unsigned int is01 = 0x01 * (1 << (dot << 1));
		unsigned int is02 = 0x02 * (1 << (dot << 1));
		unsigned int is03 = 0x03 * (1 << (dot << 1));
		unsigned int mask = 0x03 * (1 << (dot << 1));
		unsigned int nmask = ~mask;

		// Check every possible location
		for (unsigned int index = 0; index < size; index++)
		{
			unsigned int storage_i = storage_lookup[index];
			if (storage_i == -1)
				continue;

			if (dot < dots - 1)
			{
				for (unsigned int dotb = dot + 1; dotb < dots; dotb++)
				{
					// Create the bitmasks
					unsigned int is00b = 0x00 * (1 << (dotb << 1));
					unsigned int is01b = 0x01 * (1 << (dotb << 1));
					unsigned int is02b = 0x02 * (1 << (dotb << 1));
					unsigned int is03b = 0x03 * (1 << (dotb << 1));
					unsigned int maskb = 0x03 * (1 << (dotb << 1));
					unsigned int nmaskb = ~maskb;

					unsigned int otherindex;
					unsigned int storage_otheri;
					storage_location location;
					if (((index & mask) == is01) && ((index & maskb) == is02b))
					{
						// Add +tunnel_rates at(dota, dotb) where dota is 0x01 and dotb is 0x02 and turns to 0x00 and 0x03 (and opposite), others don't care
						otherindex = (index & nmask & nmaskb) | is00 | is03b;
						storage_otheri = storage_lookup[otherindex];
						location = { storage_i , storage_otheri };
						tunnel_couplings_ip[dots * dot + dotb].push_back(location);
						otherindex = (index & nmask & nmaskb) | is03 | is00b;
						storage_otheri = storage_lookup[otherindex];
						location = { storage_i , storage_otheri };
						tunnel_couplings_ip[dots * dot + dotb].push_back(location);
					}
					else if (((index & mask) == is02) && ((index & maskb) == is01b))
					{
						// Add -tunnel_rates at(dota, dotb) where dota is 0x02 and dotb is 0x01 and turns to 0x00 and 0x03 (and opposite), others don't care
						otherindex = (index & nmask & nmaskb) | is00 | is03b;
						storage_otheri = storage_lookup[otherindex];
						location = { storage_i , storage_otheri };
						tunnel_couplings_im[dots * dot + dotb].push_back(location);
						otherindex = (index & nmask & nmaskb) | is03 | is00b;
						storage_otheri = storage_lookup[otherindex];
						location = { storage_i , storage_otheri };
						tunnel_couplings_im[dots * dot + dotb].push_back(location);
					}
				}
			}

			if ((index & mask) == is01)
			{
				// Add + larmor_frequencies[dot] / 2 where this dot is 0x02, others don't care
				larmor_frequencies_im[dot].push_back(storage_i);
			}
			else if ((index & mask) == is02)
			{
				// Add - larmor_frequencies[dot] / 2 where this dot is 0x01, others don't care
				larmor_frequencies_ip[dot].push_back(storage_i);
			}

			if ((index & mask) == is03)
			{
				// Add charging_energies[dot] where this dot is 0x03, others don't care
				charging_energies_i[dot].push_back(storage_i);
			}

			if ((index & mask) == is01)
			{
				// Add microwave_control[dot] where this dot is 0x01 and turns to 0x02 (and vice versa), others don't care
				unsigned int otherindex = (index & nmask) | is02;
				unsigned int storage_otheri = storage_lookup[otherindex];
				storage_location location = { storage_i, storage_otheri };
				microwave_control_i[dot].push_back(location);
			}

			if ((index & mask) == is03)
			{
				// Add - detuning_control[dot] where this dot is 0x03, others don't care
				detuning_control_im[dot].push_back(storage_i);
				for (unsigned int d = 0; d < dots; d++)
					if (d != dot)
						detuning_control_ip[d].push_back(storage_i);
			}

			// Store locations required to determine the measurement probabilities
			if ((index & mask) == is00)
				probability_i[dot * 4 + 0].push_back(storage_i);
			else if ((index & mask) == is01) {
				unsigned int otherindex = (index & nmask) | is02;
				unsigned int storage_otheri = storage_lookup[otherindex];
				probability_i[dot * 4 + 1].push_back(storage_i);
				probability_i[dot * 4 + 2].push_back(storage_otheri);
			}
			else
				probability_i[dot * 4 + 3].push_back(storage_i);
		}
	}
}

spine::systems::system_n_spin_n_singlet::~system_n_spin_n_singlet()
{
	delete[] storage_lookup;
	delete[] probability_i;
	delete[] charging_energies_i;
	delete[] larmor_frequencies_ip;
	delete[] larmor_frequencies_im;
	delete[] tunnel_couplings_ip;
	delete[] tunnel_couplings_im;
	delete[] microwave_control_i;
	delete[] detuning_control_ip;
	delete[] detuning_control_im;
}

bool spine::systems::system_n_spin_n_singlet::isValid(unsigned int index)
{
	if (index >= size)
		return false;

	// Count the number of electrons
	unsigned int electrons = 0;
	for (unsigned int i = 0; i < dots; i++)
	{
		unsigned int code = index & 0x03;
		if ((code == 0x01) || (code == 0x02))
			electrons += 1;
		else if (code == 0x03)
			electrons += 2;
		index >>= 2;
	}

	// The number of electrons must match the number of quantum dots
	return (electrons == dots);
}

int spine::systems::system_n_spin_n_singlet::getIndex(unsigned int * states)
{
	unsigned int index = 0;
	for (unsigned int i = 0; i < dots; i++)
		index += (states[i] << (i << 1));
	return storage_lookup[index];
}

vector<unsigned int> * spine::systems::system_n_spin_n_singlet::getIndexMeasurement(unsigned int dot, unsigned int state)
{
	return &probability_i[4 * dot + state];
}

unsigned int spine::systems::system_n_spin_n_singlet::getDimension()
{
	return storage_size;
}

void spine::systems::system_n_spin_n_singlet::updateHamiltonian(realnum * H)
{
	// Clear the entire Hamiltonian
	for (unsigned int i = 0; i < (storage_size * storage_size); i++)
		H[i] = 0;

	// Add the charging energies
	for (unsigned int dot = 0; dot < dots; dot++) {
		realnum value = charging_energy[dot];
		for (vector<unsigned int>::iterator it = charging_energies_i[dot].begin(); it != charging_energies_i[dot].end(); ++it) {
			unsigned int i = *it;
			H[i * storage_size + i] += value;
		}
	}

	// Add the larmor frequencies
	for (unsigned int dot = 0; dot < dots; dot++) {
		realnum value = larmor_frequency[dot];
		for (vector<unsigned int>::iterator it = larmor_frequencies_ip[dot].begin(); it != larmor_frequencies_ip[dot].end(); ++it) {
			unsigned int i = *it;
			H[i * storage_size + i] += value / 2;
		}
		for (vector<unsigned int>::iterator it = larmor_frequencies_im[dot].begin(); it != larmor_frequencies_im[dot].end(); ++it) {
			unsigned int i = *it;
			H[i * storage_size + i] -= value / 2;
		}
	}

	// Add the tunnel couplings
	for (unsigned int dot = 0; dot < dots - 1; dot++) {
		for (unsigned int dotb = dot + 1; dotb < dots; dotb++) {
			realnum value = tunnel_control[dot * dots + dotb];
			for (vector<storage_location>::iterator it = tunnel_couplings_ip[dot * dots + dotb].begin(); it != tunnel_couplings_ip[dot * dots + dotb].end(); ++it) {
				storage_location i = *it;
				H[i.x * storage_size + i.y] += value;
				H[i.y * storage_size + i.x] += value;
			}
			for (vector<storage_location>::iterator it = tunnel_couplings_im[dot * dots + dotb].begin(); it != tunnel_couplings_im[dot * dots + dotb].end(); ++it) {
				storage_location i = *it;
				H[i.x * storage_size + i.y] -= value;
				H[i.y * storage_size + i.x] -= value;
			}
		}
	}

	// Add the microwave control
	for (unsigned int dot = 0; dot < dots; dot++) {
		realnum value = microwave_control[dot] * rabi_frequency[dot];
		for (vector<storage_location>::iterator it = microwave_control_i[dot].begin(); it != microwave_control_i[dot].end(); ++it) {
			storage_location i = *it;
			H[i.x * storage_size + i.y] += value;
			H[i.y * storage_size + i.x] += value;
		}
	}

	// Add the detuning control
	for (unsigned int dot = 0; dot < dots; dot++) {
		realnum value = detuning_control[dot];
		for (vector<unsigned int>::iterator it = detuning_control_ip[dot].begin(); it != detuning_control_ip[dot].end(); ++it) {
			unsigned int i = *it;
			H[i * storage_size + i] += value;
		}
		for (vector<unsigned int>::iterator it = detuning_control_im[dot].begin(); it != detuning_control_im[dot].end(); ++it) {
			unsigned int i = *it;
			H[i * storage_size + i] -= value;
		}
	}
}
