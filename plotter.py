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
        print("Please provide (1) a list of files to print, followed by \n\
               (2) a mode to print in: either --Einf for just that term, \n\
                   or --full for the minimum compression rate plot.")
        return 0
    if args[-1] == '--full':
        for file in args[1:-1]:
            vals = pd.read_csv(file).values
            X = vals[:, 0]
            Y = -vals[:, 1] + X + entr(X)
            plt.plot(X, Y)
        plt.plot(X, X*0 + 1)
        plt.show()
    elif args[-1] == '--Einf':
        for file in args[1:-1]:
            vals = pd.read_csv(file).values
            X = vals[:, 0]
            Y = vals[:, 1]
            plt.plot(X, Y)
        plt.show()
    else:
        print("Invalid mode")

if __name__ == '__main__':
    main()

