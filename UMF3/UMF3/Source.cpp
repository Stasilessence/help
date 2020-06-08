#define _CRT_SECURE_NO_WARNINGS
#include "Structs.h"
#include "Matrix.h"
#include "LU.h"
#include <vector>
#include "solver.h" 
#include <stdio.h>
#include <ctime>
#include "functions.h" 
#include "fragmentationGrid.h"
#define forn(i, n)  for (int i = 0; i < n; ++i)
using namespace std;
void readPoints(const char* nm, vector<point> &p)
{
	int n = 0;
	double a, b;
	FILE *fd = fopen(nm, "r");
	fscanf(fd, "%d", &n);
	for (int i = 0; i < n; ++i)
	{
		fscanf(fd, "%lf%lf", &a, &b);
		p.push_back({ a, b });
	}
}

void readElems(const char* nm, vector<elem> &p, double sig, double lam)
{
	int n = 0, l1, l2, l3, l4, ln, sn;
	FILE* fd = fopen(nm, "r");
	fscanf(fd, "%d", &n);
	for (int i = 0; i < n; ++i)
	{
		fscanf(fd, "%d%d%d%d%d%d", &ln, &sn, &l1, &l2, &l3, &l4);
		p.emplace_back();
		p[i].vertex[0] = l1 - 1;
		p[i].vertex[1] = l2 - 1;
		p[i].vertex[2] = l3 - 1;
		p[i].vertex[3] = l4 - 1;
		p[i].lambda = [&](double x, double y, double t) {return lam; };
		p[i].sigma = [&](double x, double y, double t) {return sig; };
	}
}

void readFirst(Matrix &mt, const char *nm)
{
	int n = 0, l1;
	double l2, l3;

	FILE* fd = fopen(nm, "r");
	fscanf(fd, "%d", &n);
	for (int i = 0; i < n; ++i)
	{
		first_boundary fst;
		fscanf(fd, "%d%lf%lf", &l1, &l2, &l3);
		fst.vertex = l1 - 1;
		fst.bs = l2;
		fst.bo = l3;
		mt.addFirst(fst);
	}
}
void readSecond(Matrix& mt, const char* nm)
{
	int n = 0, n2, l1;
	double l2, l3;

	FILE* fd = fopen(nm, "r");
	fscanf(fd, "%d", &n);
	for (int i = 0; i < n; ++i)
	{
		second_boundary sec;
		fscanf(fd, "%d%d", &n2, &l1);
		sec.tettaS = getTettaS(l1);
		sec.tettaO = getTettaO(l1);
		for (int j = 0; j < n2; j++)
		{
			fscanf(fd, "%d", &l1);
			sec.vertex.push_back(l1 - 1);
		}
		mt.addSecond(sec);
	}

}

double calcError(double *x, vector <point> p)
{
	double sum = 0;
	double errorSum = 0;
	for (int i = 0; i < p.size(); ++i)
	{
		sum += (p[i].x * p[i].x * p[i].x  + p[i].y) + (p[i].x - p[i].y * p[i].y * p[i].y);
		errorSum += fabs(x[2 * i] - (p[i].x * p[i].x * p[i].x + p[i].y)) + fabs(x[2 * i + 1] - (p[i].x - p[i].y * p[i].y * p[i].y));
	}
	return errorSum / fabs(sum);
}


void toProfile(Matrix& in, ProfileMatrix& mt);
void calcX(ProfileMatrix* mt, double* b);
void toLUsq(ProfileMatrix* mt);
void LU(Matrix &mt,  double* x)
{
	ProfileMatrix mt2;
	toProfile(mt, mt2);
	toLUsq(&mt2);
	for (int i = 0; i < mt.n; ++i)
		x[i] = mt.b[i];
	calcX(&mt2, x);
}

int BSG(Matrix* mt, Parameters* pr, double* x);
int LUBSG(Matrix* mt, Parameters* pr, double* x);
int LUPreconditioningProcess(Matrix* mt, Parameters* pr, double* x);
int diagPreconditioningProcess(Matrix* mt, Parameters* pr, double* x);
Matrix* read_matrix(const char* fn)
{
	int n, nl, nu;
	FILE* fd = fopen(fn, "r");
	Matrix* mt = new Matrix;
	fscanf(fd, "%d", &n);
	mt->n = n;
	mt->ia = (int*)malloc(sizeof(int) * (n + 1));
	mt->diag = (double*)malloc(sizeof(double) * n);
	mt->b = (double*)malloc(sizeof(double) * n);
	forn(i, n)
		fscanf(fd, "%lf", mt->diag + i);
	forn(i, n + 1)
	{
		fscanf(fd, "%d", mt->ia + i);
		mt->ia[i] -= 1;
	}

	int nnu = mt->ia[n];
	mt->ja = (int*)malloc(sizeof(int) * nnu);
	mt->al = (double*)malloc(sizeof(double) * nnu);
	mt->au = (double*)malloc(sizeof(double) * nnu);

	forn(i, nnu)
	{
		fscanf(fd, "%d", mt->ja + i);
		mt->ja[i] -= 1;
	}

	forn(i, nnu)
		fscanf(fd, "%lf", mt->al + i);
	forn(i, nnu)
		fscanf(fd, "%lf", mt->au + i);
	forn(i, n)
		fscanf(fd, "%lf", mt->b + i);
	return mt;
}
int main()
{
	string dir = ".";
	fragmentationGrid A;
	A.process(dir);
	FILE *file;

	const vector <double> _xi = { 1E-11 };
	const vector <double> _sigma = { 1E+02 };
	const vector <double> _lambda = { 1E+01 };
	const vector <double> _omega = { 1E+00 };



	//freopen_s(&file, "test1.csv", "w", stdout);
	for (int i4 = 0; i4 < _omega.size(); ++i4)
	{
		for (int i3 = 0; i3 < _lambda.size(); ++i3)
		{
			for (int i2 = 0; i2 < _sigma.size(); ++i2)
			{
				for (int i1 = 0; i1 < _xi.size(); ++i1)
				{
					double *x;
					Matrix mt;
					Parameters p;
					vector<elem> els;
					vector<point> asd;
					p.e = 1e-10;
					p.k = 5000;

					readPoints("points.txt", asd);
					readElems("nodex.txt", els, _sigma[i2], _lambda[i3]);
					readFirst(mt, "first.txt");
					readSecond(mt, "second.txt");

					x = new double[100000];
					fill(x, x + 100000, 0);
					mt.kxi_c = _xi[i1];
					mt.w = _omega[i4];
					mt.sigma_c = _sigma[i2];
					mt.lambda_c = _lambda[i3];
					mt.init(els, asd);
					mt.build(1);
					mt.applySecond(0);
					mt.applyFirst();

					clock_t t1 = clock();
					//Matrix *mt2 = read_matrix("tht.txt");
				//	BSG(&mt, &p, x);
					//LUBSG(&mt, &p, x);
					LOS(&mt, &p, x);
					//LU(mt, x);
					int r = clock() - t1;
					printf("%d; %.2le\n", r, calcError(x, asd));
					fprintf_s(stderr, "%d %d %d %d\n", i1, i2, i3, i4);
				}
			}
		}
	}

	return 0;
}