# Phony targets
.PHONY: all clean release debug

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall 

# Debug and Release flags
DEBUG_FLAGS = -g -O0 -DDEBUG
RELEASE_FLAGS = -O3

# Project structure
SRC_DIR = src
OBJ_DIR = obj
INCLUDE_DIR = include

# Source files and headers
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
HEADERS = $(wildcard $(INCLUDE_DIR)/*.h)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))
EXECUTABLE = kdv

# Default mode is Release
BUILD_TYPE ?= release

# if Debug
ifeq ($(BUILD_TYPE), debug)
    CXXFLAGS += $(DEBUG_FLAGS)
else
    # if Release
    CXXFLAGS += $(RELEASE_FLAGS)
endif

# Display the flags being used for debug
#$(info CXXFLAGS=$(CXXFLAGS))

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

# Compilation rule for Release mode
release:
	$(MAKE) BUILD_TYPE=release

# Compilation rule for Debug mode
debug:
	$(MAKE) BUILD_TYPE=debug


# Cleaning build files
clean:
	rm -rf $(OBJ_DIR) $(EXECUTABLE) *.txt gnuplot.log frames/
