CC=clang
MPI_BIN_DIR=/home/felipe/Coding/c/openmpi/openmpi-4.1.1/bin
MPICC=OMPI_MPICC=$(CC) $(MPI_BIN_DIR)/mpicc
SRC=$(wildcard *.c)
TARGET=$(SRC:%.c=%)
CFLAGS=
LDFLAGS=-lm

NP ?= $(shell nproc)
FILE ?=results_4.plt

VALGRIND_COMMAND:=valgrind -q --tool=memcheck\
		--track-origins=yes

.PHONY: run debug clear clean benchmark results statistics save
all: $(TARGET)

# for compiling
$(TARGET): %: %.c
	$(MPICC) $(CFLAGS) $(LDFLAGS) $< -o $@

# for running
# to change number of processes to n, run this with the NP=n argument
run: $(TARGET)
	@for target in $^; do\
		echo "running '$$target' with $(NP) processes...\n";\
		$(MPI_BIN_DIR)/mpirun --oversubscribe -np $(NP) $$target;\
	done;

plot:
	gnuplot -e "plot '$(FILE)' using 1:2 title 'quicksort' with lines,\
		'$(FILE)' using 1:3 title 'hyperquicksort' with lines;\
		pause -1 \"Hit any key to continue\""

debug: CFLAGS+=-g -O0 -DDEBUG
debug: $(TARGET)
debug:
	@for target in $^; do\
		echo "running '$$target' with $(NP) processes...\n";\
		$(MPI_BIN_DIR)/mpirun --oversubscribe -np $(NP) $(VALGRIND_COMMAND) ./$$target;\
	done;

statistics:
	@cat $(FILE) | cut -d ' ' -f 2,3 | \
		awk '{ qsum+=$$1; hsum+=$$2 } END \
		{ speedup=qsum/hsum; print "speedup: " speedup " efficiency: " speedup/np " processes: " np  }' \
		np="$$(echo "$(FILE)" | grep -o "_[[:digit:]]*." | sed 's/[_\.]//g')"

clean: clear
	-@rm -f $(TARGET)
