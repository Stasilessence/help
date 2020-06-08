#include "fragmentationGrid.h"

void fragmentationGrid::read(string dir)
{
	FILE* f;
	fopen_s(&f, (dir + "/x.txt").c_str(), "r");
	fscanf_s(f, "%d", &howx);
	xw = new double[howx];
	ihx = new int[howx - 1];
	for (int i = 0; i < howx; ++i) {
		fscanf_s(f, "%lf", &xw[i]);
	}
	fclose(f);
	fopen_s(&f, (dir + "/y.txt").c_str(), "r");
	fscanf_s(f, "%d", &howy);
	yw = new double[howy];
	ihy = new int[howy - 1];
	for (int i = 0; i < howy; ++i) {
		fscanf_s(f, "%lf", &yw[i]);
	}
	fclose(f);
	fopen_s(&f, (dir + "/hxhy.txt").c_str(), "r");
	for (int i = 0; i < howx - 1; ++i) {
		fscanf_s(f, "%d", &ihx[i]);
		hx = (xw[i + 1] - xw[i]) / ihx[i];
		for (int j = 0; j < ihx[i]; ++j)
		{
			x.push_back(xw[i] + j * hx);
		}
	}
	x.push_back(xw[howx - 1]);
	for (int i = 0; i < howy - 1; ++i) {
		fscanf_s(f, "%d", &ihy[i]);
		hy = (yw[i + 1] - yw[i]) / ihy[i];
		for (int j = 0; j < ihy[i]; j++)
		{
			y.push_back(yw[i] + j * hy);
		}
	}
	y.push_back(yw[howy - 1]);
	xsize = x.size();
	ysize = y.size();
	fclose(f);
}

double fragmentationGrid::NodeValue1(double x, double y)
{
	return x * x * x + y;
}

double fragmentationGrid::NodeValue2(double x, double y)
{
	return x - y * y * y;
}

void fragmentationGrid::process(string dir)
{
	read(dir);
	FILE* f;

	fopen_s(&f, (dir + "/points.txt").c_str(), "w");
	fprintf(f, "%d\n", ysize * xsize);
	for (int i = 0; i < ysize; i++)
		for (int j = 0; j < xsize; j++)
			fprintf(f, " %0.16lg %0.16lg\n", x[j], y[i]);
	fclose(f);
	fopen_s(&f, (dir + "/nodex.txt").c_str(), "w");

	fprintf(f, "%d\n", (xsize - 1) * (ysize - 1));
	for (int i = 0; i < (xsize - 1) * (ysize - 1); i++) {
		int k = K(i / (xsize - 1) % 2 ? (xsize - 1) * (i / (xsize - 1) + 1) - i % (xsize - 1) - 1: i) + 1;
		fprintf(f, "0 0\n");
		fprintf(f, "%d %d ", k, k + 1);
		fprintf(f, "%d %d\n\n", k + xsize, k + xsize + 1);
	}
	fclose(f);
	fopen_s(&f, (dir + "/first.txt").c_str(), "w");
	fprintf(f, "%d\n", xsize * 2 + ysize * 2 - 4);
	for (int i = 0; i < xsize; i++)
	{
		fprintf(f, "%d %.16le %.16le\n", i + 1, NodeValue1(x[i], y[0]), NodeValue2(x[i], y[0]));
	}
	for (int i = 0; i < xsize; i++)
	{
		fprintf(f, "%d %.16le %.16le\n", i + (xsize) * (ysize - 1) + 1, NodeValue1(x[i], y[ysize - 1]), NodeValue2(x[i], y[ysize - 1]));
	}
	for (int i = 1; i < ysize - 1; i++)
	{
		fprintf(f, "%d %.16le %.16le\n", i * (xsize) + 1, NodeValue1(x[0], y[i]), NodeValue2(x[0], y[i]));
	}
	for (int i = 1; i < ysize - 1; i++)
	{
		fprintf(f, "%d %.16le %.16le\n", (i + 1) * (xsize) - 1 + 1, NodeValue1(x[xsize - 1], y[i]), NodeValue2(x[xsize - 1], y[i]));
	}
	fclose(f);

	delete[]xw;
	delete[]yw;
	delete[]ihx;
	delete[]ihy;
	x.~vector();
	y.~vector();
}

int fragmentationGrid::K(int i)
{
	return i / (xsize - 1) * xsize + i % (xsize - 1);
}
