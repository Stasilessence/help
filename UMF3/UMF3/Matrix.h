#pragma once
#include <vector>
#include "Structs.h"

using namespace std;

class Matrix
{
public:
	double *diag, *al, *au, *b, w, kxi_c, sigma_c, omega_c, lambda_c;
	int n, *ia, *ja, m;
	vector<elem> els;
	vector<point> points;
	vector<first_boundary> firsts;
	vector<second_boundary> seconds;
	void buildLocal(elem el, double tc);
	void build(double tc);
	void addFirst(first_boundary);
	void addSecond(second_boundary);
	
	void applyFirst();
	void applySecond(double tc);
	void init(vector<elem> elems, vector<point> v);

};

