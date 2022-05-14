#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
using namespace std;
 
// Iterative DP function to find the number of times
// the second string occurs in the first string,
// whether continuous or discontinuous
double count(string a, string b)
{
    int m = a.length();
    int n = b.length();
 
    // Create a table to store results of sub-problems
    double lookup[m + 1][n + 1];
    memset(lookup, 0, (m+1)*(n+1)*sizeof(double));
 
    // If first string is empty
    for (int i = 0; i <= n; ++i)
        lookup[0][i] = 0;
 
    // If second string is empty
    for (int i = 0; i <= m; ++i)
        lookup[i][0] = 1;
 
    // Fill lookup[][] in bottom up manner
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            // If last characters are same, we have two
            // options -
            // 1. consider last characters of both strings
            //    in solution
            // 2. ignore last character of first string
            if (a[i - 1] == b[j - 1])
                lookup[i][j] = lookup[i - 1][j - 1] +
                               lookup[i - 1][j];
                 
            else
                // If last character are different, ignore
                // last character of first string
                lookup[i][j] = lookup[i - 1][j];
        }
    }
 
    return lookup[m][n];
}

double rand01() {
	return ((double) rand()) / RAND_MAX;
}

string BDC(string x, double d) {
	string y = "";
	for (char c : x){
		if (rand01() < 1 - d) 
			y += c;
	}
	return y;
}

string sample_ber(int n, double d) {
	string out = "";
	for (int i = 0; i < n; i++) {
		out += rand01() < d ? "1" : "0";
	}
	return out;
}

double estimate_Einf(int n, double d, int nsamples) {
	double sum_vals = 0;
	for (int i = 0; i < nsamples; i++) {
		string x = sample_ber(n, 0.5);
		string y = BDC(x, d);

		//cout << "x = " << x <<endl;
		//cout << "y = " << y <<endl;
		double c = count(x, y);
		if (c <= 0) {
		       	cout << "Error: y should appear at least once in x" << endl;
			cout << "x = " << x << endl;
			cout << "y = " << y << endl;
			cout << "count = " << c <<endl;
		}
		sum_vals += log2((double) c) / (double) n;
	}
	return sum_vals / nsamples;
}

// Driver code
int main(int argc, char* argv[])
{
	if (argc != 4) {
		cout << "Please provide a number of values of d to plot (first arg)," << endl
		     << "a value of n (second arg), and a number of iterations (third arg)" << endl;
		return 0;
	}
	int d_num = atoi(argv[1]);
	int n = atoi(argv[2]);
	int iters = atoi(argv[3]);
	double d_incr = 1 / (double) (d_num - 1);
	
	double d_vals[d_num];
	double E_inf_vals[d_num];
	int i = 0;
	for (double d = 0; i < d_num; d += d_incr, i++) {
		d_vals[i] = d;
		double estimate = estimate_Einf(n, d, iters);
		E_inf_vals[i] = estimate;
	}
	ofstream file;
	file.open("sim_output.csv");
	for (int j = 0; j < d_num; j++) {
		file << d_vals[j] << "," << E_inf_vals[j] << endl;
	}
	file.close();
	return 0;
}
