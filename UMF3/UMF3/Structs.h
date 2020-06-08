#pragma once
#include <functional>
#include <vector>
using namespace std;
typedef function<double(double, double, double)> Func;


struct elem {
	int vertex[4];
	Func lambda, sigma, khi;
};
struct first_boundary {
	int vertex;
	double bs, bo;
};
struct second_boundary {
	vector<int> vertex;
	Func tettaS, tettaO;
};
struct point {
	double x, y;
};
typedef struct node {
	double x, y, z;
	int number;
} Node;

typedef struct tetrahedron {
	Node* nodes[4];
	double lambda, gamma;
} Tetrahedron;

typedef struct triangle {
	Node* nodes[3];
} Triangle;




typedef struct parameters {
	int k;
	double e;
} Parameters;