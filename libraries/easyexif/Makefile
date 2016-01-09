CXX=g++
CXXFLAGS=-O2 -pedantic -Wall -Wextra -ansi -std=c++11

all: demo

exif.o: exif.cpp
	$(CXX) $(CXXFLAGS) -c exif.cpp

demo: exif.o demo.cpp
	$(CXX) $(CXXFLAGS) -o demo exif.o demo.cpp

clean:
	rm -f *.o demo demo.exe
