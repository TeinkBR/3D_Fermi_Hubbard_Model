#ifndef OPERATORS_H
#define OPERATORS_H

#include <complex>

const int L = 10; // Linear size of the lattice
const int N = L*L*L; // Number of lattice sites

// Define the creation and annihilation operators for each lattice site
extern std::complex<double> c[N], c_dag[N];

#endif
