# How to run the simulations, lower bound computations, and plotter

There are 3 files of interest:
- `sim.cpp`: Implements the simulation to compute the term E_inf = lim_{n \to \infty} (1/n) \E_{X, Y}[\log_2(|DD(X, Y)|)].
- `upper.cpp`: Implements the (deterministic) upper bound on the E_inf term based on the random walk analysis.
- `plotter.py`: Plots the outputs produced by the `sim` and `upper` programs.

Below are instructions to run each of the programs above:
- `sim.cpp`: First run the command `make` in the command line to compile all files. Then run `./exec_sim`. You'll get an error saying you need to pass certain command-line arguments. Passing those arguments correctly will make the program run, and produce (or re-write) the file `sim_output.csv`.
- `upper.cpp`: First run the command `make` in the command line to compile all files. Then run `./exec_upp`. You'll get an error saying you need to pass certain command-line arguments. Passing those arguments correctly will make the program run, and produce (or re-write) the file `upper_output.csv`.
- `plotter.py`: Once you've produced the data you want to plot by runnin one or both of the programs above, run `python3 plotter.py`. You'll again be prompted for command-line arguments. Once you call this again with the correct arguments, a matplotlib interactive window will pop up with the plots.

If you get stuck on anything, needless to say you should ping Fran!
