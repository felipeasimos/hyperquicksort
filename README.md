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
	speedup: 2.79116 efficiency: 0.697791 processes: 4
MPI Statistics:
	speedup: 2.46827 efficiency: 0.617068 processes: 4

# show statistics using data from 'results_mpi_2.plt' and 'results_openmp_2.plt'
make statistics NP=2
OpenMP Statistics:
	speedup: 1.98561 efficiency: 0.992807 processes: 2
MPI Statistics:
	speedup: 1.8177 efficiency: 0.908849 processes: 2

# call gnuplot to show us the results for 4 processes/threads or 2 processes/threads
make plot
make plot NP=2
```

![Plot showing quicksort curve above the hyperquicksort curve, getting farther as both of them get higher as they go right. \
The Y Axis represents time and the X Axis is the array size](./README_results_4_plot.png)
