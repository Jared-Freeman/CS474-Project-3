# Check for OS Type
UNAME_S := $(shell uname -s)

# Linux
ifeq ($(UNAME_S), Linux)
	CC=g++
# Mac
else
	CC=clang++
	LIBS= -stdlib=libc++
endif

#Compiler Flags
CXXFLAGS=-g -Wall -std=c++0x

# .o Compilation
O_FILES=main.o image.o Mask.o

# Point to includes of local directories
INCLUDES=-I../include


all: $(O_FILES)
	$(CC) $(CXXFLAGS) -o CS474 $(O_FILES) $(LIBS)
	
image.o: ../src/image.cpp
	$(CC) $(CXXFLAGS) -c ../src/image.cpp -o image.o $(INCLUDES)
	
main.o: ../src/main.cpp
	$(CC) $(CXXFLAGS) -c ../src/main.cpp -o main.o $(INCLUDES)
	
Mask.o: ../src/Mask.cpp
	$(CC) $(CXXFLAGS) -c ../src/Mask.cpp -o Mask.o $(INCLUDES)	


clean:
	-@if rm *.o CS474 2>/dev/null || true; then echo "Main Removed"; else echo "No Main"; fi
