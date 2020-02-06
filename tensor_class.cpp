/*
 * tensor_class.cpp
 *
 *  Created on: Jan 23, 2020
 *      Author: chikitkin
 */
#include <iostream>
#include <vector>
using namespace std;

class Tensor {
public:
	// Constructor
	Tensor(){
		n1 = 0;
		n2 = 0;
		n3 = 0;
	}

	Tensor(int n1_, int n2_, int n3_){
		// Create zero tensor
		n1 = n1_;
		n2 = n2_;
		n3 = n3_;
		// full tensor
		core = new double[n1 * n2 * n3];
	}
	// Destructor
	~Tensor(){
		delete [] core;
	}
	vector<int> shape(){
		auto dim = {n1, n2, n3};
		return dim;
	}
	double& At(int i1, int i2, int i3){
		// map from 1d index and 3-multi-index:
		// I = i1 * n2 * n3 + i2 * n3 + i3;
		return core[I(i1, i2, i3)];
	}


	// rounding
	void round(){
		// do nothing for full tensor
	}
	// compute sum of all elements
	double sum() {
		double s = 0.;
		for (int i = 0; i < n1 * n2 * n3; ++i){
			s += core[i];
		}
		return s;
	}

	// element-wise summation
	friend Tensor add(Tensor& t1, Tensor& t2){
		// check, that shapes are equal;
		if (t1.shape() != t2.shape()){
			cout << "Different shapes in sum!" << endl;
			exit(-1);
		}
		vector<int> dim = t1.shape();
		Tensor result(dim[0],  dim[1] , dim[2]);

		for (int i = 0; i < dim[0] * dim[1] * dim[2]; ++i){
			result.core[i] = t1.core[i] + t2.core[i];
		}
		return result;
	}

	//element-wise multiplication
	friend Tensor mult(Tensor& t1, Tensor& t2){
		// check, that shapes are equal;
		if (t1.shape() != t2.shape()){
			cout << "Different shapes in sum!" << endl;
			exit(-1);
		}
		vector<int> dim = t1.shape();
		Tensor result(dim[0],  dim[1] , dim[2]);

		for (int i = 0; i < dim[0] * dim[1] * dim[2]; ++i){
			result.core[i] = t1.core[i] * t2.core[i];
		}
		return result;
	}

	void print(){
		// print slices
		for (int i3 = 0; i3 < n3; ++i3){
			cout << "T(:, :, " + to_string(i3) + ")" << endl;
			for (int i1 = 0; i1 < n1; ++i1){
				for (int i2 = 0; i2 < n2; ++i2){
					cout << to_string(core[I(i1, i2, i3)]) + " ";
				}
				cout << endl;
			}
		}
	}

private:
	int I(int i1, int i2, int i3){
		return i1 * n2 * n3 + i2 * n3 + i3;
	}
	vector<int> multiI(int I){
		// TODO: implement
	}
	// dynamic array to store parameters
	double* core;
	// sizes along each dimension;
	int n1, n2, n3;

};


// TODO: overload operators "+, *";
Tensor operator +(Tensor& t1, Tensor& t2){
	return add(t1, t2);
}

Tensor operator *(Tensor& t1, Tensor& t2){
	return mult(t1, t2);
}

int main(){
	int n1, n2, n3;
	n1 = 3;
	n2 = 3;
	n3 = 3;
	Tensor t(n1, n2, n3);
	for (int i1 = 0; i1 < n1; ++i1){
		for (int i2 = 0; i2 < n2; ++i2){
			for (int i3 = 0; i3 < n3; ++i3){
				t.At(i1, i2, i3) = i1 + i2 + i3;
			}
		}
	}
	t.print();
	double sum = t.sum();
	cout << "t.sum() = " << sum << endl;

	Tensor t2 = t;
	Tensor t3 = t2 + t;
	t3.print();

	t3 = t2 * t;
	t3.print();
	cout << "End" << endl;
	return 0;
}


