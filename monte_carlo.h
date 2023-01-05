#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "operators.h"

// Define the Fermi-Hubbard model parameters
const double t = 1.0; // Hopping strength
const double U = 2.0; // Onsite interaction strength

// Define the temperature and chemical potential
const double T = 1.0;
const double mu = 0.0;

// Define the number of Monte Carlo steps and the random number generator
const int nsteps = 10000;
extern std::mt19937 generator;

// Perform the Monte Carlo loop
void monte_carlo_loop();

#endif
