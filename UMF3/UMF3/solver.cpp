#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structs.h"
#include "Matrix.h"
#include <stdlib.h>
#include "assert.h"
#define forn(i, n)  for (int i = 0; i < n; ++i)

void mul_matrix_vector(Matrix * mt, double *x, double *rs)
{
	int n = mt->n, *il = mt->ia, *jl = mt->ja, *iu = mt->ia, *ju = mt->ja;
	double *lells = mt->al;
	double *uells = mt->au;
	forn(i, n)
		rs[i] = x[i] * mt->diag[i];
	forn(i, n)
	{
		double sum = 0;
		int s = il[i];
		int e = il[i + 1];
		for (int j = s; j < e; ++j)
			sum += x[jl[j]] * lells[j];
		rs[i] += sum;
	}
	forn(i, n)
	{
		double sum = 0;
		int s = iu[i];
		int e = iu[i + 1];
		int j = s;
		for (; j < e; ++j)
			rs[ju[j]] += x[i] * uells[j];
	}
}
void mul_matrix_T_vector(Matrix * mt, double *x, double *rs)
{
	int n = mt->n, *il = mt->ia, *jl = mt->ja, *iu = mt->ia, *ju = mt->ja;
	double *lells = mt->au;
	double *uells = mt->al;

	forn(i, n)
		rs[i] = x[i] * mt->diag[i];
	forn(i, n)
	{
		double sum = 0;
		int s = il[i];
		int e = il[i + 1];
		for (int j = s; j < e; ++j)
			sum += x[jl[j]] * lells[j];
		rs[i] += sum;
	}
	forn(i, n)
	{
		double sum = 0;
		int s = iu[i];
		int e = iu[i + 1];
		int j = s;
		for (; j < e; ++j)
			rs[ju[j]] += x[i] * uells[j];
	}
}
void mul_matrix_diagmatrixR(Matrix * mt, double *x)
{
	int n = mt->n, *il = mt->ia, *jl = mt->ja, *iu = mt->ia, *ju = mt->ja;
	double *lells = mt->al;
	double *uells = mt->au;
	double *diag = mt->diag;
	for (int i = 0; i < n; ++i)
	{
		double sum = 0;
		int s = il[i];
		int e = il[i + 1];
		for (int j = s; j < e; ++j)
			lells[j] = x[i] * lells[j];
		diag[i] = x[i] * diag[i];
		s = iu[i];
		e = iu[i + 1];
		for (int j = s; j < e; ++j)
			uells[j] = x[ju[j]] * uells[j];
	}
}
void mul_matrix_diagmatrixL(Matrix * mt, double *x)
{
	int n = mt->n, *il = mt->ia, *jl = mt->ja, *iu = mt->ia, *ju = mt->ja;
	double *lells = mt->al;
	double *uells = mt->au;
	double *diag = mt->diag;
	for (int i = 0; i < n; ++i)
	{
		double sum = 0;
		int s = il[i];
		int e = il[i + 1];
		for (; s < e; ++s)
			lells[s] = x[jl[s]] * lells[s];
		diag[i] = x[i] * diag[i];
		s = iu[i];
		e = iu[i + 1];
		for (; s < e; ++s)
			uells[s] = x[i] * uells[s];
	}
}
double scalar(double *x1, double *x2, int n)
{
	double rs = 0;
	forn(i, n)
		rs += x1[i] * x2[i];
	return rs;
}
int diagPreconditioningProcess(Matrix * mt, Parameters * pr, double *x)
{
	double a, b;
	int n = mt->n, i;
	double *f = mt->b;
	double *r = new double[n];
	double *z = new double[n];
	double *p = new double[n];
	double *Ar = new double[n];
	double *LU = new double[n];

	forn(i, n)
		LU[i] = sqrt(1 / mt->diag[i]);
	mul_matrix_vector(mt, x, r);
	forn(i, n)
		r[i] = LU[i] * (f[i] - r[i]);

	forn(i, n)
		z[i] = LU[i] * r[i];

	mul_matrix_diagmatrixR(mt, LU);
	mul_matrix_vector(mt, z, p);
	mul_matrix_diagmatrixL(mt, LU);
	double sqff = sqrt(scalar(f, f, n));
	double sqrr = sqrt(scalar(r, r, n));
	for (i = 0; i < pr->k && sqrr / sqff > pr->e; ++i)
	{
		double pp = scalar(p, p, n);
		a = scalar(p, r, n) / pp;
		forn(i, n)
			x[i] = x[i] + a * z[i];
		forn(i, n)
			r[i] = r[i] - a * p[i];
		mul_matrix_vector(mt, r, Ar);
		b = -scalar(p, Ar, n) / pp;
		forn(i, n)
			z[i] = LU[i] * r[i] + b * z[i];
		forn(i, n)
			p[i] = Ar[i] + b * p[i];
		sqrr = sqrt(scalar(r, r, n));
		printf("\r%d %.14e", i, sqrr / sqff);
	}
	printf("\n");
	delete LU;
	delete r;
	delete z;
	delete p;
	delete Ar;
	return i;
}
void toLUsq(Matrix * mt)
{
	int n = mt->n, *il = mt->ia, *jl = mt->ja, *iu = mt->ia, *ju = mt->ja;
	double *lells = mt->al;
	double *uells = mt->au;
	double *diag = mt->diag;
	forn(i, n)
	{
		int ls = il[i];
		int le = il[i + 1];
		for (int is = ls; is < le; ++is)
		{
			double sum = lells[is];
			for (int js = ls, ujs = iu[jl[is]]; js < is && ujs < iu[jl[is] + 1];)
			{
				if (ju[ujs] < jl[js]) {
					++ujs;
					continue;
				}
				if (ju[ujs] == jl[js]) {
					sum -= lells[js] * uells[ujs];
					++ujs;
				}
				js += 1;
			}
			lells[is] = sum / diag[jl[is]];
		}
		int us = iu[i];
		int ue = iu[i + 1];
		for (int is = us; is < ue; ++is)
		{
			double sum = uells[is];
			for (int js = us, ljs = il[ju[is]]; js < is && ljs < il[ju[is] + 1]; )
			{
				if (jl[ljs] < ju[js]) {
					++ljs;
					continue;
				}
				if (jl[ljs] == ju[js]) {
					sum -= uells[js] * lells[ljs];
					++ljs;
				}
				js += 1;
			}
			uells[is] = sum / diag[ju[is]];
		}
		double sum = diag[i];
		for (; us < ue && ls < le;)
		{
			if (ju[us] < jl[ls]) {
				++us;
				continue;
			}
			if (jl[ls] == ju[us]) {
				sum -= uells[us] * lells[ls];
				++us;
			}
			//if (i == 199)
				//printf("\t%lf  %lf\n", uells[us], lells[ls]);
			ls += 1;
		}
		//printf("\t%lf %d %lf\n", sum, i, diag[i]);
		//if (sum < 0)
		//	exit(0);
		assert(sum >= 0);
		diag[i] = sqrt(sum);
	}
}

void mul_RLmatrix_vector(Matrix * mt, double *x, double *rs)
{
	int n = mt->n, *il = mt->ia, *jl = mt->ja, *iu = mt->ia, *ju = mt->ja;
	double *lells = mt->al;
	double *uells = mt->au;
	double *diag = mt->diag;
	forn(i, n)
	{
		double sum = x[i];
		int s = il[i];
		int e = il[i + 1];
		for (int j = s; j < e; ++j)
			sum -= rs[jl[j]] * lells[j];
		rs[i] = sum / diag[i];
	}
}
void mul_RLmatrix_T_vector(Matrix * mt, double *x, double *rs)
{
	int n = mt->n, *il = mt->ia, *jl = mt->ja, *iu = mt->ia, *ju = mt->ja;
	double *lells = mt->au;
	double *uells = mt->al;
	double *diag = mt->diag;
	forn(i, n)
	{
		double sum = x[i];
		int s = il[i];
		int e = il[i + 1];
		for (int j = s; j < e; ++j)
			sum -= rs[jl[j]] * lells[j];
		rs[i] = sum / diag[i];
	}
}
void mul_RUmatrix_vector(Matrix * mt, double *x, double *rs)
{
	int n = mt->n, *il = mt->ia, *jl = mt->ja, *iu = mt->ia, *ju = mt->ja;
	double *lells = mt->al;
	double *uells = mt->au;
	double *diag = mt->diag;
	forn(i, n)
		rs[i] = x[i];
	for (int i = n - 1; i >= 0; --i)
	{
		int s = iu[i];
		int e = iu[i + 1];
		int j = s;
		rs[i] /= diag[i];
		for (; j < e; ++j)
			rs[ju[j]] -= rs[i] * uells[j];
	}
}



void mul_Umatrix_vector(Matrix * mt, double *x, double *rs)
{
	int n = mt->n, *il = mt->ia, *jl = mt->ja, *iu = mt->ia, *ju = mt->ja;
	double *lells = mt->al;
	double *uells = mt->au;
	double *diag = mt->diag;
	forn(i, n)
		rs[i] = 0;
	for (int i = n - 1; i >= 0; --i)
	{
		int s = iu[i];
		int e = iu[i + 1];
		int j = s;
		rs[i] += x[i] * diag[i];
		for (; j < e; ++j)
			rs[ju[j]] += x[i] * uells[j];
	}

}

void mul_RUmatrix_T_vector(Matrix * mt, double *x, double *rs)
{
	int n = mt->n, *il = mt->ia, *jl = mt->ja, *iu = mt->ia, *ju = mt->ja;
	double *lells = mt->au;
	double *uells = mt->al;
	double *diag = mt->diag;
	forn(i, n)
		rs[i] = x[i];
	for (int i = n - 1; i >= 0; --i)
	{
		int s = iu[i];
		int e = iu[i + 1];
		int j = s;
		rs[i] /= diag[i];
		for (; j < e; ++j)
			rs[ju[j]] -= rs[i] * uells[j];
	}

}


void copyMatrix(Matrix * s, Matrix * d)
{
	int n = s->n, *il = s->ia, *jl = s->ja, *iu = s->ia, *ju = s->ja;
	double *lells = s->al;
	double *uells = s->au;
	double *diag = s->diag;
	d->n = n;
	d->ia = (int *)malloc(sizeof(int) * (n + 1));
	d->ja = (int *)malloc(sizeof(int) * (n + 1));
	d->diag = (double *)malloc(sizeof(double) * n);
	int nnl = il[n];
	int nnu = iu[n];
	d->ja = (int *)malloc(sizeof(int) * nnl);
	//d->ja = malloc(sizeof(int) * nnu);
	d->al = (double *)malloc(sizeof(double) * nnl);
	d->au = (double *)malloc(sizeof(double) * nnu);
	forn(i, n + 1)
	{
		d->ia[i] = il[i];

	}
	forn(i, n)
		d->diag[i] = diag[i];
	forn(i, nnl)
	{
		d->ja[i] = jl[i];
		d->al[i] = lells[i];
	}
	forn(i, nnu)
	{
		d->ja[i] = ju[i];
		d->au[i] = uells[i];
	}
}

int LUPreconditioningProcess(Matrix * mt, Parameters * pr, double *x)
{
	double a, b;
	double *f = mt->b;
	int n = mt->n, i;
	Matrix *LU = (Matrix *)malloc(sizeof(Matrix));
	double *r = (double *)malloc(sizeof(double) * n);
	double *z = (double *)malloc(sizeof(double) * n);
	double *p = (double *)malloc(sizeof(double) * n);
	double *Ar = (double *)malloc(sizeof(double) * n);
	double *bf = (double *)malloc(sizeof(double) * n);
	copyMatrix(mt, LU);
	toLUsq(LU);
	mul_matrix_vector(mt, x, r);
	forn(i, n)
		r[i] = f[i] - r[i];
	mul_RLmatrix_vector(LU, r, r);

	mul_RUmatrix_vector(LU, r, z);
	mul_matrix_vector(mt, z, bf);
	mul_RLmatrix_vector(LU, bf, p);

	double sqff = sqrt(scalar(f, f, n));
	double sqrr = sqrt(scalar(r, r, n));
	for (i = 0; i < pr->k && sqrr / sqff > pr->e; ++i)
	{
		double pp = scalar(p, p, n);
		a = scalar(p, r, n) / pp;
		forn(i, n)
			x[i] = x[i] + a * z[i];
		forn(i, n)
			r[i] = r[i] - a * p[i];

		mul_RUmatrix_vector(LU, r, bf);
		mul_matrix_vector(mt, bf, Ar);
		mul_RLmatrix_vector(LU, Ar, Ar);

		b = -scalar(p, Ar, n) / pp;
		forn(i, n)
			z[i] = bf[i] + b * z[i];
		forn(i, n)
			p[i] = Ar[i] + b * p[i];
		sqrr = sqrt(scalar(r, r, n));
		printf("tap:%d %.14e\n", i, sqrr / sqff);
	}
	printf("\n");
	free(LU->al);
	free(LU->diag);
	free(LU->ia);
	free(LU->ja);
	free(LU);
	free(r);
	free(z);
	free(p);
	free(Ar);
	free(bf);
	return 0;
}


int LOS(Matrix * mt, Parameters * pr, double *x)
{
	double a, b;
	int n = mt->n, i;
	double *f = mt->b;
	double *r = new double[n];
	double *z = new double[n];
	double *p = new double[n];
	double *Ar = new double[n];

	memset(r, 0, n);
	memset(z, 0, n);
	memset(p, 0, n);
	memset(Ar, 0, n);

	mul_matrix_vector(mt, x, r);
	forn(i, n)
		r[i] = f[i] - r[i];

	forn(i, n)
		z[i] = r[i];


	mul_matrix_vector(mt, z, p);

	double sqff = sqrt(scalar(f, f, n));
	double sqrr = sqrt(scalar(r, r, n));
	for (i = 0; i < pr->k && sqrr / sqff > pr->e; ++i)
	{
		double pp = scalar(p, p, n);
		a = scalar(p, r, n) / pp;
		forn(i, n)
			x[i] = x[i] + a * z[i];
		forn(i, n)
			r[i] = r[i] - a * p[i];
		mul_matrix_vector(mt, r, Ar);
		b = -scalar(p, Ar, n) / pp;
		forn(i, n)
			z[i] = r[i] + b * z[i];
		forn(i, n)
			p[i] = Ar[i] + b * p[i];
		sqrr = sqrt(scalar(r, r, n));

	}

	printf("%d; %.12e; ", i - 1, sqrr / sqff);
	delete r;
	delete z;
	delete p;
	delete Ar;
	return i;
}


int BSG(Matrix * mt, Parameters * pr, double *x)
{
	double a, b;
	int n = mt->n, i;
	double *f = mt->b;
	double *r = new double[n];
	double *z = new double[n];
	double *p = new double[n];
	double *s = new double[n];

	double *Az = new double[n];



	mul_matrix_vector(mt, x, r);
	forn(i, n)
		r[i] = f[i] - r[i];

	forn(i, n)
		z[i] = p[i] = s[i] = r[i];





	double sqff = sqrt(scalar(f, f, n));
	double sqrr = sqrt(scalar(r, r, n));
	for (i = 0; i < pr->k && sqrr / sqff > pr->e; ++i)
	{
		mul_matrix_vector(mt, z, Az);
		double sAz = scalar(s, Az, n);

		a = scalar(p, r, n) / sAz;

		double prkm1 = scalar(p, r, n);
		forn(i, n)
			x[i] = x[i] + a * z[i];
		forn(i, n)
			r[i] = r[i] - a * Az[i];

		mul_matrix_T_vector(mt, s, Az);

		forn(i, n)
			p[i] = p[i] - a * Az[i];


		b = scalar(p, r, n) / prkm1;
		forn(i, n)
			z[i] = r[i] + b * z[i];

		forn(i, n)
			s[i] = p[i] + b * s[i];
		sqrr = sqrt(scalar(r, r, n));
		printf("tap:%d %.14e\n", i, sqrr / sqff);
	}
	printf("\n");
	delete r;
	delete z;
	delete p;
	delete Az;
	return i;
}




int LUBSG(Matrix * mt, Parameters * pr, double *x)
{
	double a, b;
	double *f = mt->b;
	int n = mt->n, i;
	Matrix *LU = (Matrix *)malloc(sizeof(Matrix));
	double *r = new double[n];
	double *z = new double[n];
	double *p = new double[n];
	double *s = new double[n];
	double *Az = new double[n];
	double *bf = new double[n];
	double *xl = new double[n];
	//double* LU = new double[n];
	copyMatrix(mt, LU);
	toLUsq(LU);

	mul_RUmatrix_vector(LU, x, r);
	mul_matrix_vector(mt, r, bf);

	forn(i, n)
		r[i] = f[i] - bf[i];

	mul_RLmatrix_vector(LU, r, s);

	forn(i, n)
		z[i] = p[i] = r[i] = s[i];

	//mul_matrix_diagmatrixR(mt, LU);
	//mul_matrix_diagmatrixL(mt, LU);

	double sqff = sqrt(scalar(f, f, n));
	double sqrr = sqrt(scalar(r, r, n));

	for (i = 0; i < pr->k && sqrr / sqff > pr->e; ++i)
	{

		mul_RUmatrix_vector(LU, z, bf);
		mul_matrix_vector(mt, bf, Az);
		mul_RLmatrix_vector(LU, Az, Az);

		double sAz = scalar(s, Az, n);
		double prkm1 = scalar(p, r, n);

		a = prkm1 / sAz;

		forn(i, n)
			x[i] = x[i] + a * z[i];

		forn(i, n)
			r[i] = r[i] - a * Az[i];


		mul_RUmatrix_T_vector(LU, s, bf);
		mul_matrix_T_vector(mt, bf, Az);
		mul_RLmatrix_T_vector(LU, Az, Az);


		forn(i, n)
			p[i] = p[i] - a * Az[i];


		b = scalar(p, r, n) / prkm1;

		forn(i, n)
			z[i] = r[i] + b * z[i];
		forn(i, n)
			s[i] = p[i] + b * s[i];


		sqrr = sqrt(scalar(r, r, n));

		
	}
	printf("%d; %.14e; ", i - 1, sqrr / sqff);
	mul_RUmatrix_vector(LU, x, x);
	/*forn(i, n)
		x[i] = x[i] * LU[i];*/


	free(LU);
	free(r);
	free(z);
	free(p);
	free(bf);
	return 0;
}

int DIAGBSG(Matrix * mt, Parameters * pr, double *x)
{
	double a, b;
	double *f = mt->b;
	int n = mt->n, i;
	//	Matrix* LU = (Matrix*)malloc(sizeof(Matrix));
	double *r = new double[n];
	double *z = new double[n];
	double *p = new double[n];
	double *s = new double[n];
	double *Az = new double[n];
	double *bf = new double[n];
	double *xl = new double[n];
	double *LU = new double[n];
	//	copyMatrix(mt, LU);
	//	toLUsq(LU);
	forn(i, n)
		LU[i] = sqrt(1 / mt->diag[i]);

	forn(i, n)
		x[i] = x[i] * LU[i];
	//	mul_RUmatrix_vector(LU, x, x);


	mul_matrix_vector(mt, x, r);

	forn(i, n)
		r[i] = LU[i] * (f[i] - r[i]);

	//mul_RLmatrix_vector(LU, r, r);


	forn(i, n)
		z[i] = p[i] = s[i] = r[i];

	mul_matrix_diagmatrixR(mt, LU);
	mul_matrix_diagmatrixL(mt, LU);

	double sqff = sqrt(scalar(f, f, n));
	double sqrr = sqrt(scalar(r, r, n));
	for (i = 0; i < pr->k && sqrr / sqff > pr->e; ++i)
	{

		mul_matrix_vector(mt, z, Az);


		double sAz = scalar(s, Az, n);
		double prkm1 = scalar(p, r, n);
		a = prkm1 / sAz;

		forn(i, n)
			x[i] = x[i] + a * z[i];

		forn(i, n)
			r[i] = r[i] - a * Az[i];




		mul_matrix_T_vector(mt, s, Az);


		forn(i, n)
			p[i] = p[i] - a * Az[i];


		b = scalar(p, r, n) / prkm1;

		forn(i, n)
			z[i] = r[i] + b * z[i];
		forn(i, n)
			s[i] = p[i] + b * s[i];


		sqrr = sqrt(scalar(r, r, n));

		printf("tap:%d %.14e\n", i, sqrr / sqff);
	}
	//mul_RUmatrix_vector(LU, x, x);
	forn(i, n)
		x[i] = x[i] * LU[i];
	printf("\n");


	free(LU);
	free(r);
	free(z);
	free(p);
	free(bf);
	return 0;
}