# Phony targets
.PHONY: all clean

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2

# Project structure
SRC_DIR = src
OBJ_DIR = obj
INCLUDE_DIR = include

# Source files and headers
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
HEADERS = $(wildcard $(INCLUDE_DIR)/*.h)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))
EXECUTABLE = kdv

# Default target
all: $(EXECUTABLE)

# Linking the executable
$(EXECUTABLE): $(OBJECTS) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compiling object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -c $< -o $@

# Directories
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR) 


# Cleaning build files
clean:
	rm -rf $(OBJ_DIR) $(EXECUTABLE)
