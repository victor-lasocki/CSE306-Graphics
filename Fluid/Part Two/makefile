CXX = clang++
CXXFLAGS = -g3 -Wall -Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include
LDFLAGS = -L/usr/local/opt/libomp/lib -lomp -llbfgs

all: fluid

clean:
	rm -f fluid main.o

fluid: main.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^  # Include LDFLAGS for linking

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
