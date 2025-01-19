.PHONY: all clean

#Compiler
CC = g++
CFLAGS = -Wall 

#Executable name
TARGET = proj

#Source files
SRC = proj.cpp

#Compilation rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

#Clean up
clean:
	rm -f $(TARGET) 
