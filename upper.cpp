#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;

double tget(vector<vector<double>>& table, int i, int j) {
	return table[table.size() / 2 + i][j];
}

void tset(vector<vector<double>>& table, int i, int j, double val) {
	table[table.size() / 2 + i][j] = val;
}

void update(vector<vector<double>>& table, int n, int k, double d) {
	double p1 = tget(table, k - 1, n - 1);
	double p0 = tget(table, k, n - 1);
	double p_1 = tget(table, k + 1, n - 1);
	if (k != 0)
		tset(table, k, n, (1 - d) * pow(2, -3.0 / 2) * p1 + 0.5 * ((1-d)/2 + d)* p0 + d * pow(2, -3.0/2) * p_1);
	else
		tset(table, k, n, (1 - d) * pow(2, -3.0/2) * p1 + p0 / 2 + d * pow(2, -3.0/2) * p_1);

}
double get_nth_upper(int N, double d) {
	vector<vector<double>> table;
	for (int i = 0; i < 2 * N + 3; i++) {
		vector<double> row;
		for (int j = 0; j < N; j++)
			row.push_back(0);
		table.push_back(row); 
	}
	tset(table, 0, 0, 1);
	for (int n = 1; n < N; n++) {
		for (int k = -N; k <= N; k++) {
			update(table, n, k, d);
		}
	}
	double res = tget(table, 0, N - 1);
	// std::cout << res << endl;
	return 1 + log2(res) / N;
}


// Driver code
int main(int argc, char* argv[])
{
	if (argc != 3) {
		cout << "Please provide a number of values of d to compute (first arg), and" << endl
		     << "a value of n (second arg)." << endl;
		return 0;
	}
	int d_num = atoi(argv[1]);
	int n = atoi(argv[2]);
	double d_incr = 1 / (double) (d_num - 1);

	double d_vals[d_num];
	double E_inf_vals[d_num];

	int i = 0;
	for (double d = 0; i < d_num; d += d_incr, i++) {
		d_vals[i] = d;
		double estimate = get_nth_upper(n, d);
		E_inf_vals[i] = estimate;
	}
	ofstream file;
	file.open("upper_output.csv");
	for (int j = 0; j < d_num; j++) {
		file << d_vals[j] << "," << E_inf_vals[j] << endl;
	}
	file.close();
	return 0;
}
