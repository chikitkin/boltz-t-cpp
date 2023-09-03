#ifndef HEADER_H_
#define HEADER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <memory>

#include <omp.h>

#include "mkl.h"
#include "mkl_lapacke.h"

typedef double REAL;

// void MKL_imatcopy(char ordering, char trans, size_t rows, size_t cols, REAL *alpha, REAL *a, size_t src_lda, size_t dst_lda) {
//     return MKL_Dimatcopy(ordering, trans, rows, size_t cols, *alpha, *a, src_lda, dst_lda)
// }

#endif /* HEADER_H_ */
