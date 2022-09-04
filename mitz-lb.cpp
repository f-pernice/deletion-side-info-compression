#include <iostream>
#include <string.h>
#include <assert.h>
#include <cmath>
#include <vector>
#include <fstream>
using namespace std;

double entr(double p) {
	const double eps = 1e-20;
	if (abs(p) < eps || abs(1 - p) < eps) return 0;
	return -p * log2(p) - (1 - p) * log2(1 - p);
}

long double nChoosek( long double n, long double k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    long double result = n;
    for( long double i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}
long double get_term(int k, int i, int z, int r, int s, long double d) {
	long double term1 = powl((1-d)/d, k);
	long double term2 = nChoosek(z+r, k) - nChoosek(r, k);
	long double term3 = powl(d / 2, z + r + s);
	long double term4;
	if (i == 0) {
		if (r > 0 || s > 0) term4 = 0;
		else term4 = 1;
	} else term4 = nChoosek(r-1, i-1) * nChoosek(s-1, i-1);

	if (term2 < 1) cout << "term less than 1!!" << endl;
	long double term5 = log2l(term2);
	if (term1 < 0) cout << "negative term 1" << endl;
	if (term2 < 0) cout << "negative term 2" << endl;
	if (term3 < 0) cout << "negative term 3" << endl;
	if (term4 < 0) cout << "negative term 4" << endl;
	if (term5 < 0) cout << "negative term 5" << endl;
	
	return ((1-d)/2) * term1*term2*term3*term4*term5;
}


int main(int argc, char* argv[]) {
	if (argc != 7) {
		cout << "You provided " << argc << " arguments." << endl;
		cout << "Exactly 6 arguments needed:" << endl;
		cout << "(1) a value of kmax >= 1" << endl;
		cout << "(2) a value of imax >= 0" << endl;
		cout << "(3) a value of zmax >= 1" << endl;
		cout << "(4) a value of rmax >= imax" << endl;
		cout << "(5) a value of smax >= imax" << endl;
		cout << "(6) a value of dvals >= 1" << endl;
		return 0;
	}
	int kmax = atoi(argv[1]);
	int imax = atoi(argv[2]);
	int zmax = atoi(argv[3]);
	int rmax = atoi(argv[4]);
	int smax = atoi(argv[5]);
	int dvals = atoi(argv[6]);
	double d_incr = 1 / (double) (dvals - 1);

	assert(kmax >= 1);
	assert(imax >= 0);
	assert(zmax >= 1);
	assert(rmax >= imax);
	assert(smax >= imax);
	assert(dvals >= 1);
	vector<long double> d_values;
	vector<long double> Y_values;

	int diter = 0;
	for (long double d = 0; diter < dvals; d += d_incr, diter++) {
		cout << "On iter " << diter << " out of " << dvals << endl;
		long double sum = 0;
		d_values.push_back(d);
		for (int i = 0; i <= imax; i++) {
			for (int z = 1; z <= zmax; z++) {
				for (int r = i; r <= rmax + i; r++) {
					for (int s = i; s <= smax + i; s++) {
						for (int k = 1; k <= min(kmax, r + z); k++) {
							sum += get_term(k, i, z, r, s, d);
						}
					}
				}
			}
		}
		Y_values.push_back(sum);
	}

	ofstream f;
	f.open("mitz-lb_" + to_string(kmax) + "_" + to_string(imax) + "_" + to_string(zmax) + "_" + to_string(rmax) + "_" + to_string(smax) + "_" + to_string(dvals) + ".csv");
	for (int i = 0; i < (int) d_values.size(); i++) {
		f << d_values[i] << "," << Y_values[i] << endl;
	}
	f.close();
	return 0;
}
