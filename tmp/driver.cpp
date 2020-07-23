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
using namespace std;

void driver(const string& path, bool write_restart = 1) {

	ofstream logfile;
	logfile.open ("log.txt");
	logfile << path << endl;
	logfile.close();

	if (write_restart) {


	}

}

