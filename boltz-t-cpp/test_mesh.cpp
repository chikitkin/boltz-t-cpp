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
using namespace std;

int main() {

	vector < vector < int > > A;

	vector <int> a;
	a.push_back(1);
	a.push_back(2);
	a.push_back(3);
	a.push_back(4);
	a.push_back(5);
	a.push_back(6);

	int x1 = a[0];
	int x2 = a[1];
	int x3 = a[2];
	int x4 = a[3];
	int x5 = a[4];
	int x6 = a[5];

	a[3] = x1;
	a[4] = x6;
	a[5] = x4;
	a[0] = x2;
	a[1] = x5;
	a[2] = x3;

	for (int n : a) {
		cout << n << endl;
	}


	cout << "a.capacity " << a.capacity() << endl;

	A.push_back(vector <int> {1, 2, 3, 4, 5});
	A.push_back(vector <int> {1, 2, 3, 4, 5});
	A.push_back(vector <int> {1, 2, 3, 4, 5});
	A.push_back(a);
	A.push_back(vector <int> {1, 2, 3, 4, 5});

	cout << "A[3].capacity " << A[3].capacity() << endl;

	A.shrink_to_fit();

	return 0;
}
