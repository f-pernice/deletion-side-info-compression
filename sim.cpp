#include <iostream>
#include <string.h>
#include <assert.h>
#include <cmath>
#include <fstream>
using namespace std;
#define EPS 1e-10
 
// Iterative DP function to find the number of times
// the second string occurs in the first string,
// whether continuous or discontinuous
long double count(string a, string b)
{
    int m = a.length();
    int n = b.length();
 
    // Create a table to store results of sub-problems
    long double **lookup = new long double*[m + 1];
    for (int i = 0; i <= m; i++) {
    	lookup[i] = new long double[n + 1];
	for (int j = 0; j <= n; j++)
		lookup[i][j] = 0;
    }
 
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
    long double answer = lookup[m][n];
    for (int i = 0; i <= m; i++) delete[] lookup[i];
    delete[] lookup; 
    cout << "count = " << answer << endl;
    return answer; 
}

long double rand01() {
	return ((long double) rand()) / RAND_MAX;
}

string BDC(string x, long double d) {
	string y = "";
	for (char c : x){
		if (rand01() < 1 - d) 
			y += c;
	}
	return y;
}
/*
string BDC_fixed_string_start(string x, long double d) {
	if (abs(d) < EPS) return x;
	if (abs(1-d) < EPS) return "";
	string y = "";
	int incr = ceil(1/(1-d));
//	cout << "incr = " << incr << endl;
	int dels = 0;
	for (int i = 0; i < x.length() && dels < d * x.length(); i++) {
		if (i % incr == 0) y += x[i];
		else dels++;
	}
//	cout << "y = " + y << endl;
	return y;
}*/

string sample_ber(int n, long double d) {
	string out = "";
	for (int i = 0; i < n; i++) {
		out += rand01() < d ? "1" : "0";
	}
	return out;
}

long double estimate_Einf(int n, long double d, int nsamples, bool uncorrelated_y) {
	long double sum_vals = 0;
	for (int i = 0; i < nsamples; i++) {
		string x = sample_ber(n, 0.5);
		string y;
		if (!uncorrelated_y) y = BDC(x, d);
//		if (!uncorrelated_y) y = BDC_fixed_string_start(x, d);
		else y = sample_ber((int) ((1 - d) * n), 0.5);

		long double c = count(x, y);
		if (!uncorrelated_y && c <= 0) {
		       	cout << "Error: y should appear at least once in x" << endl;
			cout << "x = " << x << endl;
			cout << "y = " << y << endl;
			cout << "count = " << c <<endl;
		}
		if (!uncorrelated_y) sum_vals += log2l(c) / (long double) n;
		else sum_vals += log2l(c + 1) / (long double) n;

	}
	return sum_vals / nsamples;
}

long double estimate_Einf(int n, long double d, int nsamples) {
	return estimate_Einf(n, d, nsamples, false);
}

// Driver code
int main(int argc, char* argv[])
{
	if (!(argc == 4 || argc == 5)) {
		cout << "Please provide a number of values of d to plot (first arg)," << endl
		     << "a value of n (second arg), and a number of iterations (third arg)." << endl
		     << "You can also optionally specify a --uncorrelated flag to sample X and Y independently" << endl;
		return 0;
	}
	int d_num = atoi(argv[1]);
	int n = atoi(argv[2]);
	int iters = atoi(argv[3]);
	bool uncorrelated = false;
	if (argc == 5) {
		assert(strcmp(argv[4], "--uncorrelated") == 0);
		uncorrelated = true;
	}
	long double d_incr = 1 / (long double) (d_num - 1);
	
	long double d_vals[d_num];
	long double E_inf_vals[d_num];
	int i = 0;
	for (long double d = 0; i < d_num; d += d_incr, i++) {
		if (d <= 0.5) continue;
		d_vals[i] = d;
		cout << "On d-value number " << i << " / " << d_num << "." << endl;
		long double estimate = estimate_Einf(n, d, iters, uncorrelated);
		E_inf_vals[i] = estimate;
		cout << "Finished d-value " << i << ". Output was " << E_inf_vals[i] << "." << endl;
	}
	ofstream file;
	if (!uncorrelated) file.open("sim_output_n" + to_string(n) + ".csv");
	else file.open("sim_output_uncorrelated_n" + to_string(n) + ".csv");

	for (int j = 0; j < d_num; j++) {
		file << d_vals[j] << "," << E_inf_vals[j] << endl;
	}
	file.close();
	return 0;
}
