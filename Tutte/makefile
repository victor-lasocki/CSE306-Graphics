# Define the compiler
CXX = g++

# Define compiler flags
CXXFLAGS = -Wall -Wextra -std=c++11

# Define the target executable
TARGET = tutte

# Define the source files
SRC = main.cpp

# Define the object files
OBJ = $(SRC:.cpp=.o)

# Rule to build the target
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

# Rule to build the object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: clean
