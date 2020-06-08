#include "Matrix.h"
#include "Matrixes.h"
#include <unordered_map>
#include <algorithm>



double localMatrix[64];
double bl[8];
void Matrix::addFirst(first_boundary b)
{
	firsts.push_back(b);
}
void Matrix::addSecond(second_boundary b)
{
	seconds.push_back(b);
}
void Matrix::applySecond(double tc)
{
	for (int i = 0; i < seconds.size(); ++i)
	{
		for (int j = 1; j < seconds[i].vertex.size(); ++j)
		{
			int p1 = seconds[i].vertex[j - 1];
			int p2 = seconds[i].vertex[j];
			point v1 = points[p1];
			point v2 = points[p2];
			double tetta1S = seconds[i].tettaS(v1.x, v1.y, tc);
			double tetta1O = seconds[i].tettaO(v1.x, v1.y, tc);
			double tetta2S = seconds[i].tettaS(v2.x, v2.y, tc);
			double tetta2O = seconds[i].tettaO(v2.x, v2.y, tc);
			double h = (v1.x - v2.x != 0) ? abs(v1.x - v2.x) : abs(v1.y - v2.y);
			p1 *= 2;
			p2 *= 2;
			b[p1] += h * (2 * tetta1S + tetta2S) / 6;
			b[p1 + 1] += h * (2 * tetta1O + tetta2O) / 6;
			b[p2] += h * (tetta1S + 2 * tetta2S) / 6;
			b[p2 + 1] += h * (tetta1O + 2 * tetta2O) / 6;
		}
	}
}
void Matrix::applyFirst()
{
	for (int i = 0; i < firsts.size(); ++i)
	{
		first_boundary th = firsts[i];
		int ps = th.vertex * 2;
		b[ps] = th.bs;
		b[ps + 1] = th.bo;
		diag[ps] = 1;
		diag[ps + 1] = 1;
		int s = ia[ps];
		int e = ia[ps + 2];
		for (; s < e; ++s)
			al[s] = 0;
		for (int i = 0; i < m; ++i)
		{
			if (ja[i] == ps || ja[i] == ps + 1)
				au[i] = 0;
		}
		
	}
}
void Matrix::build(double tc)
{
	for (int i = 0; i < els.size(); ++i)
	{
		els[i].khi = [&](double x, double y, double t) {return kxi_c;};
		buildLocal(els[i], tc);
		for (int j = 0; j < 4; ++j)
		{
			int cth = els[i].vertex[j] * 2;
			for (int jj = 0; jj < 4; jj++)
			{
				int jth = els[i].vertex[jj] * 2;
				if (jth <= cth)
					continue;
				int is = 0;
				int s = ia[jth];
				int e = ia[jth + 1];
				for (; s < e && ja[s] < cth; ++s, ++is)
					;
				if (ja[s] != cth)
					continue;

				al[s] += localMatrix[jj * 16 + j * 2];
				al[s + 1] += localMatrix[jj * 16 + j * 2 + 1];
				al[e + is] += localMatrix[(jj*2 + 1) * 8 + j * 2];
				al[e + is + 1] += localMatrix[(jj*2 + 1) * 8 + j * 2 + 1];

				au[s] += localMatrix[j * 16 + jj * 2];
				au[e + is] += localMatrix[j * 16 + jj * 2 + 1];
				au[s + 1] += localMatrix[(j*2 + 1) * 8 + jj * 2];
				au[e + is + 1] += localMatrix[(j*2 + 1) * 8 + jj * 2 + 1];
			}
			
			diag[cth] += localMatrix[j * 16 + j * 2];
			diag[cth + 1] += localMatrix[j * 16 + j * 2];
			int e = ia[cth + 2] - 1;
			al[e] += -localMatrix[j * 16 + j * 2 + 1];
			au[e] += localMatrix[j * 16 + j * 2 + 1];
			b[cth] += bl[j * 2];
			b[cth + 1] += bl[j * 2 + 1];
		}

	}
}
void Matrix::buildLocal(elem el, double tc)
{
	double l1 = el.lambda(points[el.vertex[0]].x, points[el.vertex[0]].y, tc);
	double l2 = el.lambda(points[el.vertex[1]].x, points[el.vertex[1]].y, tc);
	double l3 = el.lambda(points[el.vertex[2]].x, points[el.vertex[2]].y, tc);
	double l4 = el.lambda(points[el.vertex[3]].x, points[el.vertex[3]].y, tc);
	double khi1 = el.khi(points[el.vertex[0]].x, points[el.vertex[0]].y, tc);
	double khi2 = el.khi(points[el.vertex[1]].x, points[el.vertex[1]].y, tc);
	double khi3 = el.khi(points[el.vertex[2]].x, points[el.vertex[2]].y, tc);
	double khi4 = el.khi(points[el.vertex[3]].x, points[el.vertex[3]].y, tc);
	double sigma1 = el.sigma(points[el.vertex[0]].x, points[el.vertex[0]].y, tc);
	double sigma2 = el.sigma(points[el.vertex[1]].x, points[el.vertex[1]].y, tc);
	double sigma3 = el.sigma(points[el.vertex[2]].x, points[el.vertex[2]].y, tc);
	double sigma4 = el.sigma(points[el.vertex[3]].x, points[el.vertex[3]].y, tc);

	Func Fs = [&](double x, double y, double t) {
		return -6 * lambda_c * x - omega_c * sigma_c * (x - y * y * y) - omega_c * omega_c * kxi_c * (x * x * x + y);
	};
	Func Fo = [&](double x, double y, double t) {
		return 6 * lambda_c * y + omega_c * sigma_c * (x * x * x + y) - omega_c * omega_c * kxi_c * (x - y * y * y);
	};

	double fs1 = Fs(points[el.vertex[0]].x, points[el.vertex[0]].y, tc);
	double fs2 = Fs(points[el.vertex[1]].x, points[el.vertex[1]].y, tc);
	double fs3 = Fs(points[el.vertex[2]].x, points[el.vertex[2]].y, tc);
	double fs4 = Fs(points[el.vertex[3]].x, points[el.vertex[3]].y, tc);
	double fo1 = Fo(points[el.vertex[0]].x, points[el.vertex[0]].y, tc);
	double fo2 = Fo(points[el.vertex[1]].x, points[el.vertex[1]].y, tc);
	double fo3 = Fo(points[el.vertex[2]].x, points[el.vertex[2]].y, tc);
	double fo4 = Fo(points[el.vertex[3]].x, points[el.vertex[3]].y, tc);
	/////////////////////////////////////
	double hx = abs(points[el.vertex[0]].x - points[el.vertex[1]].x);
	double hy = abs(points[el.vertex[0]].y - points[el.vertex[2]].y);
	////////////////////////////////////
	for (int i = 0; i < 4; i += 1)
	{
		int li = i * 2;
		for (int j = 0; j < 4; ++j)
		{
			int lj = j * 2;
			int mtpos44 = i * 4 + j;
			localMatrix[li * 8 + lj] = localMatrix[(li + 1) * 8 + lj + 1] =
				hy / hx * (l1 * w12dx[mtpos44] + l2 * w12dx[mtpos44] + l3 * w34dx[mtpos44] + l4 * w34dx[mtpos44])
				+
				hx / hy * (l1 * w13dy[mtpos44] + l2 * w24dy[mtpos44] + l3 * w13dy[mtpos44] + l4 * w24dy[mtpos44])
				-
				w * w * hx * hy * (khi1 * w1ww[mtpos44] + khi2 * w2ww[mtpos44] + khi3 * w3ww[mtpos44] + khi4 * w4ww[mtpos44]);
			localMatrix[li * 8 + lj + 1] = localMatrix[(li + 1) * 8 + lj] = 
				w * hx * hy * (sigma1 * w1ww[mtpos44] + sigma2 * w2ww[mtpos44] + sigma3 * w3ww[mtpos44] + sigma4 * w4ww[mtpos44]);
			localMatrix[li * 8 + lj + 1] = -localMatrix[li * 8 + lj + 1];
		}
		bl[li] = hx * hy / 36 * (fs1 * C[i * 4] + fs2 * C[i * 4 + 1] + fs3 * C[i * 4 + 2] + fs4 * C[i * 4 + 3]);
		bl[li + 1] = hx * hy / 36 * (fo1 * C[i * 4] + fo2 * C[i * 4 + 1] + fo3 * C[i * 4 + 2] + fo4 * C[i * 4 + 3]);
	}
}
void tryAddVertex(vector <vector <int>> &s, int a, int b)
{
	a = a * 2;
	b = b * 2;
	auto g = find(s[a].begin(), s[a].end(), b);
	if (g != s[a].end())
		return;
	s[a].push_back(b);
	s[a].push_back(b + 1);
	s[a + 1].push_back(b);
	s[a + 1].push_back(b + 1);
}

void Matrix::init(vector<elem> elems, vector<point> v)
{
	els = elems;
	points = v;
	n = 2 * (v.size());

	diag = new double[n];
	b = new double[n];

	vector <vector <int>> set_(n);
	vector <int> rsJ, rsI;
	for (int i = 0; i < elems.size(); ++i)
	{
		int *vtx = elems[i].vertex;

		for (int ii = 0; ii < 4; ii++)
		{
			int cth = vtx[ii];
			if (set_[cth*2 + 1].size() == 0)
				set_[cth*2 + 1].push_back(cth*2);
			for (int j = 0; j < 4; ++j)
			{
				if (vtx[j] < cth)
					tryAddVertex(set_, cth, vtx[j]);
			}
		}
	}
	rsI.push_back(0);
	for (int i = 0, c = 0; i < n; ++i)
	{
		sort(set_[i].begin(), set_[i].end());
		rsJ.insert(rsJ.end(), set_[i].begin(), set_[i].end());
		c += set_[i].size();
		rsI.push_back(c);
	}
	m = rsJ.size();
	ja = new int[m];
	al = new double[m];
	au = new double[m];
	fill(al, al + m, 0);
	fill(diag, diag + n, 0);
	fill(b, b + n, 0);
	fill(au, au + m, 0);
	ia = new int[rsI.size()];
	for (int i = 0; i < n + 1; ++i)
		ia[i] = rsI[i];
	for (int i = 0; i < m; ++i)
		ja[i] = rsJ[i];
	return;
}