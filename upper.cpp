#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;

long double tget(long double **table, int i, int j, int len) {
	return table[len/ 2 + i][j];
}

void tset(long double **table, int i, int j, long double val, int len) {
	table[len / 2 + i][j] = val;
}

long double tget3D(long double ***table, int k, int j, int u, int len) {
	return table[len/ 2 + k][j][u];
}

void tset3D(long double ***table, int k, int j, int u, long double val, int len) {
	table[len/ 2 + k][j][u] = val;
}

void update(long double **table, int n, int k, long double d, bool sym, int len) {
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

void update3D(long double ***table, int j, int k, int u, long double d, int space_range_len) {
	long double D11 = tget3D(table, k, j-1, u-1, space_range_len);
	long double D10 = tget3D(table, k-1, j-1, u-1, space_range_len);
	long double D01 = tget3D(table, k+1, j-1, u, space_range_len);
	long double D00 = tget3D(table, k, j-1, u, space_range_len);
	
	tset3D(table, k, j, u, 
		d*d * D11 + d * (1-d) * (1/sqrt(2)) * D10 + d * (1-d) * (1/sqrt(2)) * D01 + (1-d)*(1-d)*(k == 0 ? 1 : 0.5) * D00,
		space_range_len); 
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

	if (!sym) return 1 + log2(res) / (N-1);
	else return entr(d) + log2(res) / (N - 1);
}

double get_nth_upper_3D(int N, double d) {
	int space_range_len = 2 * N + 3;
	int time_range_len = N;
	int weight_range_len = N;
	long double ***table = new long double **[space_range_len];
	for (int i = 0; i < space_range_len; i++) {
		long double **time_frame = new long double *[time_range_len];
		for (int j = 0; j < time_range_len; j++) {
			long double *weight_frame = new long double[weight_range_len];
			for (int k = 0; k < weight_range_len; k++)
				weight_frame[k] = 0;
			time_frame[j] = weight_frame;
		}
		table[i] = time_frame;
	}
	tset3D(table, 0, 0, 0, 1, space_range_len);

	for (int j = 1; j < N; j++) { // iterate forward in time
		for (int k = -N; k <= N; k++) { // range over space
			for (int u = 0; u < weight_range_len; u++) // range over weight
				update3D(table, j, k, u, d, space_range_len);
		}
	}
	//cout << "d = " << d << ", n = " << N-1 << ", dn = " << (int) d*(N-1) << endl;
	long double res = tget3D(table, 0, N - 1, (int) (d * (N-1)), space_range_len);
	
	for (int i = 0; i < space_range_len; i++) {
		for (int j = 0; j < time_range_len; j++) {
			delete[] table[i][j];
		}
		delete[] table[i];
	}
	delete[] table;
	return entr(d) + log2(res) / (N-1);
}

// Driver code
int main(int argc, char* argv[])
{
	if (!(argc == 3 || argc == 4)) {
		cout << "Please provide a number of values of d to compute (first arg), and" << endl
		     << "a value of n (second arg). You can also optionally specify a --sym flag to" << endl
		     << "calculate the symmetric random walk upper bound, or a --3D flag to calculate the 3D" << endl
		     << "upper bound" <<endl;
		return 0;
	}
	int d_num = atoi(argv[1]);
	int n = atoi(argv[2]);
	double d_incr = 1 / (double) (d_num - 1);
	bool sym = false; 
	bool is_3D = false;
	if (argc == 4 && strcmp(argv[3], "--3D") == 0) is_3D = true;
	if (argc == 4 && strcmp(argv[3], "--sym") == 0) sym = true;
	assert(argc == 3 || (strcmp(argv[3], "--sym") == 0 || strcmp(argv[3], "--3D") == 0));

	double d_vals[d_num];
	double E_inf_vals[d_num];

	int i = 0;
	for (double d = 0; i < d_num; d += d_incr, i++) {
		d_vals[i] = d;
		if (is_3D) {
			cout << "Computing 3D upper bound..." << endl;
			E_inf_vals[i] = get_nth_upper_3D(n, d);

		}
		else E_inf_vals[i] = get_nth_upper(n, d, sym);
	}
	ofstream file;
	if (is_3D) {
		file.open("upper_output_3D_n" + to_string(n) + ".csv");
	} else {
		if (!sym) file.open("upper_output_n" + to_string(n) + ".csv");
		else  file.open("upper_output_symmetric_n" + to_string(n) + ".csv");
	}

	for (int j = 0; j < d_num; j++) {
		file << d_vals[j] << "," << E_inf_vals[j] << endl;
	}
	file.close();
	return 0;
}
