CC=clang
MPI_BIN_DIR=/home/felipe/Coding/c/openmpi/openmpi-4.1.1/bin
MPICC=OMPI_MPICC=$(CC) $(MPI_BIN_DIR)/mpicc
SRC_MPI=$(wildcard *_mpi.c)
SRC_OPENMP=$(wildcard *_openmp.c)
TARGET_MPI=$(SRC_MPI:%_mpi.c=%_mpi)
TARGET_OPENMP=$(SRC_OPENMP:%.c=%_openmp)
TARGET=$(TARGET_MPI) $(TARGET_OPENMP)
CFLAGS=
LDFLAGS=-lm

NP ?= $(shell nproc)
FILE ?=results_4.plt

VALGRIND_COMMAND:=valgrind -q --tool=memcheck\
		--track-origins=yes

.PHONY: run debug clean statistics plot
all: $(TARGET)

# for compiling
$(TARGET_MPI): %: %.c
	$(MPICC) $(CFLAGS) $(LDFLAGS) $< -o $@

$(TARGET_OPENMP): CFLAGS += -fopenmp
$(TARGET_OPENMP): %: %.c
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

# for running
# to change number of processes to n, run this with the NP=n argument
mpi_run: $(TARGET_MPI)
	@for target in $^; do\
		echo "running '$$target' with $(NP) processes...\n";\
		$(MPI_BIN_DIR)/mpirun --oversubscribe -np $(NP) $$target;\
	done;

plot:
	gnuplot -e "plot '$(FILE)' using 1:2 title 'quicksort' with lines,\
		'$(FILE)' using 1:3 title 'hyperquicksort MPI' with lines;\
		pause -1 \"Hit any key to continue\""

mpi_debug: CFLAGS+=-g -O0 -DDEBUG
mpi_debug: $(TARGET_MPI)
mpi_debug:
	@for target in $^; do\
		echo "running '$$target' with $(NP) processes...\n";\
		$(MPI_BIN_DIR)/mpirun --oversubscribe -np $(NP) $(VALGRIND_COMMAND) ./$$target;\
	done;

statistics:
	@cat $(FILE) | cut -d ' ' -f 2,3 | \
		awk '{ qsum+=$$1; hsum+=$$2 } END \
		{ speedup=qsum/hsum; print "speedup: " speedup " efficiency: " speedup/np " processes: " np  }' \
		np="$$(echo "$(FILE)" | grep -o "_[[:digit:]]*." | sed 's/[_\.]//g')"

clean:
	-@rm -f $(TARGET)
