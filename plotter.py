import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
mpl.use('tkagg')
import pandas as pd
import sys

EPS = 1e-10

def entr(X):
    zeros_or_ones = (X < EPS) + (X > 1 - EPS)
    # Make all zeros or ones into 1/2, and then subtract that off at the end. 
    # This prevents errors that are raised when calculating log of 0.
    X = X * (1 - zeros_or_ones) + zeros_or_ones / 2
    res = -X * np.log2(X) - (1-X) * np.log2(1-X)
    return res - zeros_or_ones

def main():
    args = sys.argv

    if len(args) < 3:
        print("Please provide (1) a list of files to print, or the keyword \n\
              \"best-bounds\" for a selection of pre-recorded bounds, followed by \n\
               (2) a mode to print in: either --Einf for just that term, \n\
                   or --full for the minimum compression rate plot.")
        return 0
    if len(args) == 3 and args[1] == "best-bounds":
        selection = ["my_heuristic_lower.csv", "my_heuristic_upper.csv", 
                "sim_output_n10000.csv", "upper_output.csv", "upper_output_symmetric_n2000.csv",
                "upper_output_3D_n1000.csv"]
#        selection = ["sim_output_n10000_low_iters.csv", "upper_output.csv"]
        args = [args[0]] + selection + [args[-1]]
    if args[-1] == '--full':
        for file in args[1:-1]:
            vals = pd.read_csv(file).values
            X = vals[:, 0]
            Y = -vals[:, 1] + X + entr(X)
            plt.plot(X, Y, label=file[:file.find('.')])
        plt.plot(X, entr(X * (X < 0.5) + 1/2 * (X >= 0.5)), label="h(d) upper bound")
        plt.plot(X, X*0 + 1)
    elif args[-1] == '--Einf':
        for file in args[1:-1]:
            vals = pd.read_csv(file).values
            X = vals[:, 0]
            Y = vals[:, 1]
#            plt.plot(X, Y)
            plt.plot(X, Y, label=file[:file.find('.')])
        plt.plot(X, X * (X <= 1/2) + (entr(X) + X - 1) * (X > 1/2), label="h(d) lower bound")
    else:
        print("Invalid mode")
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()

