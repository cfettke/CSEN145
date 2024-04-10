CC = gcc
CFLAGS = -O3 -fopenmp
LIBS = -lm
DEPS =
OBJ = matmult_omp.o


%.o: %.c $(DEPS)    
	$(CC) $(CFLAGS) -c -o $@ $<


run.out: $(OBJ)
	$(CC) $(CFLAGS) $(LIBS) -o $@ $^


