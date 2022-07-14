#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;

long double tget(long double **table, int i, int j, int len) {
	return table[len/ 2 + i][j];
}

void tset(long double **table, int i, int j, double val, int len) {
	table[len / 2 + i][j] = val;
}

void update(long double **table, int n, int k, double d, bool sym, int len) {
	long double p1 = tget(table, k - 1, n - 1, len);
	long double p0 = tget(table, k, n - 1, len);
	long double p_1 = tget(table, k + 1, n - 1, len);
	if (!sym) {
		if (k != 0)
			tset(table, k, n, 
				(1 - d) * pow(2, -3.0 / 2) * p1 + 0.5 * ((1-d)/2 + d)* p0 + d * pow(2, -3.0/2) * p_1,
				len);
		else
			tset(table, k, n, 
				(1 - d) * pow(2, -3.0/2) * p1 + p0 / 2 + d * pow(2, -3.0/2) * p_1, 
				len);
	} else {
		if (k != 0)
			tset(table, k, n, 
				(d*(1-d) / pow(2, 0.5)) * (p_1 + p1) + (pow(d, 2) + pow(1 - d, 2) / 2.0) * p0,
				len);
		else
			tset(table, k, n, 
				(d*(1-d) / pow(2, 0.5)) * (p_1 + p1) + (pow(d, 2) + pow(1-d, 2)) * p0,
				len);
	}
}

double entr(double p) {
	const double eps = 1e-20;
	if (abs(p) < eps || abs(1 - p) < eps) return 0;
	return -p * log2(p) - (1 - p) * log2(1 - p);
}

double get_nth_upper(int N, double d, bool sym) {
	int len = 2 * N + 3;
	long double **table = new long double *[len];
	for (int i = 0; i < len; i++) {
		long double *row = new long double [N];
		for (int j = 0; j < N; j++)
			row[j] = 0;
		table[i] = row; 
	}
	tset(table, 0, 0, 1, len);
	for (int n = 1; n < N; n++) {
		for (int k = -N; k <= N; k++) {
			update(table, n, k, d, sym, len);
		}
	}
	long double res = tget(table, 0, N - 1, len);

	for (int i = 0; i < 2 * N + 3; i++)
		delete[] table[i];
	delete[] table;

	if (!sym) return 1 + log2(res) / N;
	else return entr(d) + log2(res) / N;
}


// Driver code
int main(int argc, char* argv[])
{
	if (!(argc == 3 || argc == 4)) {
		cout << "Please provide a number of values of d to compute (first arg), and" << endl
		     << "a value of n (second arg). You can also optionally specify a --sym flag to" << endl
		     << "calculate the symmetric random walk upper bound." << endl;
		return 0;
	}
	int d_num = atoi(argv[1]);
	int n = atoi(argv[2]);
	double d_incr = 1 / (double) (d_num - 1);
	bool sym = argc == 4;
	assert(argc == 3 || strcmp(argv[3], "--sym") == 0);

	double d_vals[d_num];
	double E_inf_vals[d_num];

	int i = 0;
	for (double d = 0; i < d_num; d += d_incr, i++) {
		d_vals[i] = d;
		E_inf_vals[i] = get_nth_upper(n, d, sym);
	}
	ofstream file;
	if (!sym) file.open("upper_output.csv");
	else  file.open("upper_output_symmetric.csv");
	for (int j = 0; j < d_num; j++) {
		file << d_vals[j] << "," << E_inf_vals[j] << endl;
	}
	file.close();
	return 0;
}
