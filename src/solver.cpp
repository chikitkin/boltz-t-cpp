#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "math.h"
#include <algorithm>
#include <numeric>
#include "tensor_class.h"
#include "mesh.h"
#include "solver.h"
using namespace std;

GasParams::GasParams() {
	Na = 6.02214129e+23; // Avogadro constant
	kB = 1.381e-23; // Boltzmann constant, J / K
	Ru = 8.3144598; // Universal gas constant

	Mol = 40e-3;
	Rg = Ru / Mol; // = self.Ru  / self.Mol  # J / (kg * K)
	m = Mol / Na; // # kg

	Pr = 2.0 / 3.0;

	C = 144.4;
	T_0 = 273.11;
	mu_0 = 2.125e-05;

	g = 5. / 3.; // # specific heat ratio
	d = 3418e-13; // # diameter of molecule
}

double GasParams::mu_suth(double T) {
	return mu_0 * ((T_0 + C) / (T + C)) * ((T / T_0) ** (3.0 / 2.0));
}

double GasParams::mu(double T) {
	return mu_suth(200.0) * ((T / 200.0) ** (0.734));
}
