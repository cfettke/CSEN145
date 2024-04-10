CC = g++
CFLAGS = -O3 -fopenmp
LIBS = -lm 
DEPS =
OBJ = generateTree.o


%.o: %.c $(DEPS)    
	$(CC) $(CFLAGS) -c -o $@ $<


run.out: $(OBJ)
	$(CC) $(CFLAGS) $(LIBS) -o $@ $^


