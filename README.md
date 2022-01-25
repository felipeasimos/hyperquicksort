# Hyperquicksort

## Requirements

* MPI
* OpenMP (comes with clang)
* Makefile
* Gnuplot (to see the results in a chart)

## Run

```
# show debug information, also detect if final array is actually sorted and detect any memory errors using valgrind (MPI itself shows some)
make clean debug
make clean run
```

You can change the number of processes/threads with the `NP` option (default is 4):

```
# change number of processes/threads using the NP variable (default is 4)
make mpi_run NP=8
make omp_run NP=2
```

Results will be written to a `results_mpi_NP.plt` or `results_openmp_NP.plt` file, where `NP` is the
number of processes/threads.

## See Results

All data is queried from the `.plt` files, generated automatically when you execute the
program. The format is simply columns of numbers, where the first is the array size,
the second is the quicksort time and the last is the hyperquicksort time.

Statistics are calculated using the quicksort time calculated in `hyperquicksort_openmp.c`,
and the hyperquisort times measured in `hyperquicksort_openmp.c` and `hyperquicksort_mpi.c`.

You can change which set of files the results will be based on with the `NP` option (it picks
files according to the number of processes/threads used to generate them).

```
# show statistics using data from 'results_mpi_4.plt' and 'results_openmp_4.plt'
make statistics
OpenMP Statistics:
	speedup: 2.52637 efficiency: 0.631591 processes: 4
	at 500000:
		speedup: 2.6081
		efficiency: 0.652026
	at 1000000:
		speedup: 2.87959
		efficiency: 0.719896
MPI Statistics:
	speedup: 2.23171 efficiency: 0.557927 processes: 4
	at 500000:
		speedup: 2.07705
		efficiency: 0.519262
	at 1000000:
		speedup: 2.4268
		efficiency: 0.6067
# show statistics using data from 'results_mpi_2.plt' and 'results_openmp_2.plt'
make statistics NP=2
OpenMP Statistics:
	speedup: 1.90582 efficiency: 0.952909 processes: 2
	at 500000:
		speedup: 2.05442
		efficiency: 1.02721
	at 1000000:
		speedup: 1.70876
		efficiency: 0.854381
MPI Statistics:
	speedup: 1.90337 efficiency: 0.951686 processes: 2
	at 500000:
		speedup: 1.89938
		efficiency: 0.949688
	at 1000000:
		speedup: 1.95389
		efficiency: 0.976947
# call gnuplot to show us the results for 4 processes/threads or 2 processes/threads
make plot
make plot NP=2
```

![Plot showing quicksort curve above the hyperquicksort curve, getting farther as both of them get higher as they go right. \
The Y Axis represents time and the X Axis is the array size](./README_results_4_plot.png)
