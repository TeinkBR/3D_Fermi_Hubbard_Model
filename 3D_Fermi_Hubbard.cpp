#include <iostream>
#include <random>
#include <complex>

const int L = 10; // Linear size of the lattice
const int N = L*L*L; // Number of lattice sites

// Define the Fermi-Hubbard model parameters
const double t = 1.0; // Hopping strength
const double U = 2.0; // Onsite interaction strength

// Define the temperature and chemical potential
const double T = 1.0;
const double mu = 0.0;

// Define the number of Monte Carlo steps and the random number generator
const int nsteps = 10000;
std::mt19937 generator(0);

// Define the creation and annihilation operators for each lattice site
std::complex<double> c[N], c_dag[N];

int main() {
    // Initialize the lattice and the auxiliary fields
    int lattice[N] = {0}; // Occupation number of each lattice site
    double delta[N] = {0}; // Auxiliary field

    // Initialize the creation and annihilation operators
    for (int i = 0; i < N; i++) {
        c[i] = std::complex<double>(0,0);
        c_dag[i] = std::complex<double>(0,0);
    }

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
            delta[i] = -delta[i];
            energy += delta_e;
            double_occupancy += 2 * lattice[i] - 1;
            magnetization += 2 * lattice[i] - 1;

            // Update the creation and annihilation operators
            if (lattice[i] == 0) {
                c[i] = std::complex<double>(1,0);
                c_dag[i] = std::complex<double>(0,0);
            } else {
                c[i] = std::complex<double>(0,0);
                c_dag[i] = std::complex<double>(1,0);
            }
        }
    }

