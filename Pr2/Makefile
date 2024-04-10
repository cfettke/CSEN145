CXX = g++
CXXFLAGS = -O3 -fopenmp
LIBS = -lm
DEPS =
OBJ =accumulator.o


%.o: %.cpp $(DEPS)    
	$(CXX) $(CXXFLAGS) -c -o $@ $<


run.out: $(OBJ)
	$(CXX) $(CXXFLAGS) $(LIBS) -o $@ $^
