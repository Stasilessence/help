#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structs.h"
#include "Matrix.h"

void mul_matrix_vector(Matrix* mt, double* x, double* rs);
void mul_matrix_diagmatrixR(Matrix* mt, double* x);
void mul_matrix_diagmatrixL(Matrix* mt, double* x);
double scalar(double* x1, double* x2, int n);
int diagPreconditioningProcess(Matrix* mt, Parameters* pr, double* x);
int LOS(Matrix* mt, Parameters* pr, double* x);
int LUPreconditioningProcess(Matrix* mt, Parameters* pr, double* x);

int LUBSG(Matrix *mt, Parameters *pr, double *x);
