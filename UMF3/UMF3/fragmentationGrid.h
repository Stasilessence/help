#pragma once

#include <vector>
#include <string>
using namespace std;



class fragmentationGrid
{
	int *ihx, *ihy;
	int howx;
	int howy;
	int xsize;
	int ysize;
	double hx; double hy;
	double* xw;
	double* yw;
	vector<double> x;
	vector<double> y;
	void read(string dir);
	double NodeValue1(double x, double y);
	double NodeValue2(double x, double y);
	int K(int i);
public:
	void process(string dir);
};
