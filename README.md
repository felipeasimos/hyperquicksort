# Hyperquicksort

## Requirements

* MPI
* Makefile
* Gnuplot (to see the results in a chart)

## Run

```
# show debug information, also detect if final array is actually sorted and detect any memory errors using valgrind (MPI itself shows some)
make clean debug
make clean run
```

You can change the number of processes with the `NP` option (default is 4):

```
# change number of processes using the NP variable (default is 4)
make run NP=8
make run NP=2
```

Results will be written to a `results_NP.plt` file, where `NP` is the
number of processes.

## See Results

All data is queried from the `.plt` files, generated automatically when you execute the
program. The format is simply columns of numbers, where the first is the array size,
the second is the quicksort time and the last is the hyperquicksort time.

```
# show statistics using data from 'results_4.plt'
make statistics
speedup: 3.23323 efficiency: 0.808309 processes: 4

# call gnuplot to show us the results
make plot
```

You can change the file where data comes from using the `FILE` option (default is `results_4.plt`):

```
# show statistics using the 'results_2.plt' file
make statistics FILE=results_2.plt
speedup: 1.97168 efficiency: 0.985842 processes: 2
```

![Plot showing quicksort curve above the hyperquicksort curve, getting farther as both of them get higher as they go right. \
The Y Axis represents time and the X Axis is the array size](./README_results_4.png)