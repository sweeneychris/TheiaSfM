all: demo

exif.o: exif.cpp
	g++ -O2 -pedantic -Wall -Wextra -ansi -c exif.cpp

demo: exif.o demo.cpp
	g++ -O2 -Wall -pedantic -Wextra -ansi -o demo exif.o demo.cpp

clean:
	rm -f *.o demo
