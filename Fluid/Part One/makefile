# Compiler and flags
CXX = /usr/local/opt/llvm/bin/clang++
CXXFLAGS = -g3 -Wall -Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include -I/usr/local/include
LDFLAGS = -L/usr/local/opt/libomp/lib -L/usr/local/lib -lomp -llbfgs

# Target executable
TARGET = voronoi

# Source files
SRCS = main.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Build the target
all: $(TARGET)

# Link the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(TARGET) $(OBJS)
