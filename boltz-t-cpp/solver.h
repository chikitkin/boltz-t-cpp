#ifndef SOLVER_H_
#define SOLVER_H_

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
using namespace std;

class GasParams {
public:
    double Na; // Avogadro constant
    double kB; // Boltzmann constant, J / K
    double Ru; // Universal gas constant

    double Mol; // = Mol
    double Rg; // = self.Ru  / self.Mol  # J / (kg * K)
    double m; // # kg

    double Pr;

    double C;
	double T_0;
	double mu_0;

	double mu_suth(double T);

	double mu(double T);

	double g; // # specific heat ratio
	double d; // # diameter of molecule
};

class Solution {
private:


public:


};

#endif /* SOLVER_H_ */
