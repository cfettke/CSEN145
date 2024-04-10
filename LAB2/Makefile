CC = g++
CFLAGS = -g -std=c++11
LIBS = -lm
DEPS = 
OBJ = gdb_tut.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

run.out: $(OBJ)
	$(CC) $(CFLAGS) $(LIBS) -o $@ $^

clean:
	rm *.o
	rm *.out
