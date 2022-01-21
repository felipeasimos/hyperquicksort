CC=clang
MPI_BIN_DIR=/home/felipe/Coding/c/openmpi/openmpi-4.1.1/bin
MPICC=OMPI_MPICC=$(CC) $(MPI_BIN_DIR)/mpicc
SRC_MPI=$(wildcard *_mpi.c)
SRC_OMP=$(wildcard *_openmp.c)
TARGET_MPI=$(SRC_MPI:%_mpi.c=%_mpi)
TARGET_OMP=$(SRC_OMP:%_openmp.c=%_openmp)
TARGET=$(TARGET_MPI) $(TARGET_OMP)
CFLAGS=
LDFLAGS=-lm

NP ?= $(shell nproc)
FILE_MPI =results_mpi_$(NP).plt
FILE_OMP =results_openmp_$(NP).plt

VALGRIND_COMMAND:=valgrind -q --tool=memcheck\
		--track-origins=yes

.PHONY: run debug clean statistics plot
all: $(TARGET)

# for compiling
$(TARGET_MPI): %: %.c
	$(MPICC) $(CFLAGS) $(LDFLAGS) $< -o $@

$(TARGET_OMP): CFLAGS += -fopenmp
$(TARGET_OMP): %: %.c
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

# for running
# to change number of processes to n, run this with the NP=n argument
mpi_run: $(TARGET_MPI)
	@for target in $^; do\
		echo "running '$$target' with $(NP) processes...\n";\
		$(MPI_BIN_DIR)/mpirun --oversubscribe -np $(NP) $$target;\
	done;

# for running
# to change number of processes to n, run this with the NP=n argument
omp_run: CFLAGS+=-fopenmp
omp_run: $(TARGET_OMP)
	@for target in $<; do\
		echo "running '$$target' with $(NP) processes...\n";\
		OMP_NUM_THREADS=$(NP) ./$$target;\
	done;

plot:
	gnuplot -e "plot '$(FILE_OMP)' using 1:2 title 'quicksort' with lines lc rgb 'blue',\
		'$(FILE_OMP)' using 1:3 title 'hyperquicksort OMP' with lines lc rgb 'red',\
		'$(FILE_MPI)' using 1:3 title 'hyperquicksort MPI' with lines lc rgb 'green';\
		pause -1 \"Hit any key to continue\""

mpi_debug: CFLAGS+=-g -O0 -DDEBUG
mpi_debug: $(TARGET_MPI)
	@for target in $^; do\
		echo "running '$$target' with $(NP) processes...\n";\
		$(MPI_BIN_DIR)/mpirun --oversubscribe -np $(NP) $(VALGRIND_COMMAND) ./$$target;\
	done;

omp_debug: CFLAGS+=-fopenmp -g -O0 -DDEBUG
omp_debug: $(TARGET_OMP)
	@for target in $^; do\
		echo "running '$$target' with $(NP) processes...\n";\
		OMP_NUM_THREADS=$(NP) $(VALGRIND_COMMAND) ./$$target;\
	done;

statistics:
	@cat $(FILE_OMP) | cut -d ' ' -f 2,3 | \
		awk '{ qsum+=$$1; hsum+=$$2 } END \
		{ speedup=qsum/hsum; print "OpenMP Statistics:\n\tspeedup: " speedup " efficiency: " speedup/np " processes: " np }' \
		np="$$(echo "$(FILE_OMP)" | grep -o "openmp_[[:digit:]]*." | sed 's/[openmp_\.]//g')"
	@cat $(FILE_MPI) | cut -d ' ' -f 2,3 | \
		awk '{ hsum+=$$2 } END \
		{ speedup=qsum/hsum; print "MPI Statistics:\n\tspeedup: " speedup " efficiency: " speedup/np " processes: " np }' \
		np="$$(echo "$(FILE_MPI)" | grep -o "mpi_[[:digit:]]*." | sed 's/[mpi_\.]//g')" \
		qsum="$$(cat $(FILE_OMP) | cut -d ' ' -f 2,3 | awk '{ qsum+=$$1 } END { print qsum }')"

clean:
	-@rm -f $(TARGET)
