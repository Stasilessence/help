#include "LU.h"
#include "Matrix.h"
#include <stdlib.h>
#include <math.h>
#include <list>
using namespace std;
int setIndexedItem(ProfileMatrix* mt, int i, int j, double it)
{
	int* low = mt->indl, * up;
	up = (mt->indu == NULL) ? mt->indl : mt->indu;
	if (i == j)
	{
		mt->diag[i] = it;
		return 0;
	}
	else if (i > j&& i - low[i + 1] + low[i] <= j)
	{
		mt->al[low[i + 1] - i + j] = it;
		return 0;
	}
	else if (j > i&& j - up[j + 1] + up[j] <= i)
	{
		mt->au[up[j + 1] - j + i] = it;
		return 0;
	}
	return 1;
}
double getItem(ProfileMatrix* mt, int i, int j)
{
	int* low = mt->indl, * up;
	up = (mt->indu == NULL) ? mt->indl : mt->indu;
	if (i == j)
		return mt->diag[i];
	else if (i > j&& i - low[i + 1] + low[i] <= j)
		return mt->al[low[i + 1] - i + j];
	else if (j > i&& j - up[j + 1] + up[j] <= i)
		return mt->au[up[j + 1] - j + i];
	return 0;
}
void toFullMatrix(ProfileMatrix* mt, double* out)
{
	for (int i = 0; i < mt->n; ++i)
		for (int j = 0; j < mt->n; ++j)
			out[i * mt->n + j] = getItem(mt, i, j);
}
void toLUsq(ProfileMatrix* mt)
{
	int* up = (mt->indu == NULL) ? mt->indl : mt->indu;
	int* lw = mt->indl;
	double* au = mt->au;
	double* al = mt->al;
	double* aij = NULL;
	for (int i = 0; i < mt->n; ++i)
	{
		double sum;
		int lj = lw[i + 1];
		int ll = lj - lw[i];
		int ls = i - lj + lw[i];
		aij = al + lw[i + 1] - i + ls;
		for (int j = ls; j < i; ++j, ++aij)
		{
			sum = aij[0];
			int uj = up[j + 1];
			int us = j - uj + up[j];
			int k = us > ls ? us : ls;
			for (int ui = uj - j + k, li = lj - i + k; ui < uj; ++ui, ++li)
				sum -= au[ui] * al[li];
			aij[0] = sum / mt->diag[j];
		}
		int uj = up[i + 1];
		int us = i - uj + up[i];
		int ul = uj - up[i];
		aij = au + up[i + 1] - i + us;
		for (int j = us; j < i; ++j, ++aij)
		{
			sum = aij[0];
			int lj = lw[j + 1];
			int ls = j - lj + lw[j];
			int k = us > ls ? us : ls;
			for (int ui = uj - i + k, li = lj - j + k; li < lj; ++li, ++ui)
				sum -= au[ui] * al[li];
			aij[0] = sum / mt->diag[j];
		}
		ll = (ul > ll) ? ll : ul;
		sum = mt->diag[i];
		for (int li = lj - ll, ui = uj - ll; li < lj; ++li, ++ui)
			sum -= al[li] * au[ui];
		mt->diag[i] = sqrt(sum);
	}
}

void toProfile(Matrix &in, ProfileMatrix &mt)
{
	mt.n = in.n;
	mt.diag = new double[in.n];
	mt.indl = mt.indu = new int[in.n + 1];
	list<double> l, u, ia;
	ia.push_back(0);
	for (int i = 0; i < in.n; ++i)
	{
		mt.diag[i] = in.diag[i];
		for (int j = in.ia[i]; j < in.ia[i + 1] - 1; ++j)
		{
			int ls = in.ja[j];
			int rs = in.ja[j + 1];
			l.push_back(in.al[j]);
			u.push_back(in.au[j]);
			for (int k = 0; k < rs - ls - 1; ++k)
			{
				l.push_back(0);
				u.push_back(0);
			}
		}
		if (in.ia[i] < in.ia[i + 1])
		{
			l.push_back(in.al[in.ia[i + 1] - 1]);
			u.push_back(in.au[in.ia[i + 1] - 1]);
			for (int k = 0; k < i - in.ja[in.ia[i + 1] - 1] - 1; ++k)
			{
				l.push_back(0);
				u.push_back(0);
			}
		}
		ia.push_back(l.size());
	}
	mt.al = new double[l.size()];
	mt.al = (double *)memset(mt.al, 0, sizeof(double) * l.size());
	mt.au = new double[u.size()];
	mt.au = (double *)memset(mt.au, 0, sizeof(double) * u.size());
	for (int i = 0; !ia.empty(); ++i)
	{
		mt.indl[i] = ia.front();
		ia.pop_front();
	}
	for (int i = 0; !l.empty(); ++i)
	{
		mt.al[i] = l.front();
		mt.au[i] = u.front();
		l.pop_front();
		u.pop_front();
	}
	return;
}

void calcX(ProfileMatrix* mt, double* b)
{
	int* up = (mt->indu == NULL) ? mt->indl : mt->indu;
	int* lw = mt->indl;
	double* au = mt->au;
	double* al = mt->al;
	for (int i = 0; i < mt->n; ++i)
	{
		int lj = lw[i + 1];
		int ll = i - lj + lw[i];
		double sum = b[i];
		for (int j = ll, li = lj - i + ll; j < i; ++j, ++li)
			sum -= (double)al[li] * b[j];
		b[i] = sum / mt->diag[i];
	}
	for (int i = mt->n - 1; i >= 0; --i)
	{
		int uj = up[i + 1];
		int ul = i - uj + up[i];
		b[i] /= mt->diag[i];
		for (int j = ul, ui = uj - i + ul; j < i; ++j, ++ui)
			b[j] -= au[ui] * b[i];

	}
}