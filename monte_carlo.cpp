#include "monte_carlo.h"
#include <random>

std::mt19937 generator(0);

void monte_carlo_loop() {
    // Initialize the lattice and the auxiliary fields
    int lattice[N] = {0}; // Occupation number of each lattice site
    double delta[N] = {0}; // Auxiliary field

    // Initialize the observables
    double energy = 0;
    double double_occupancy = 0;
    double magnetization = 0;

    // Perform the Monte Carlo loop
    for (int step = 0; step < nsteps; step++) {
        // Select a lattice site randomly
        std::uniform_int_distribution<int> distribution(0, N-1);
        int i = distribution(generator);

        // Calculate the energy change due to flipping the occupation of the site
        double delta_e = U * lattice[i] * (lattice[i] - 1) - 2 * t * delta[i] * (lattice[i] - 0.5);

        // Calculate the acceptance probability
        double p = std::exp(-delta_e / T);

        // Perform the Monte Carlo step
        if (p > 1 or p > distribution(generator)) {
            // Accept the move
            lattice[i] = 1 - lattice[i];
            delta[i] =
