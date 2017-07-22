CXXFLAGS = -std=c++11 -Wall -O2 -Wno-parentheses -Wno-misleading-indentation
CPPFLAGS = -I../include
LDFLAGS =
LIBS =

HEADERS = $(wildcard ../include/MatOp/*.h) $(wildcard ../include/LinAlg/*.h) $(wildcard ../include/Util/*.h) $(wildcard ../include/*.h)


.PHONY: all test clean

all: QR.out Eigen.out SymEigs.out SymEigsShift.out GenEigs.out GenEigsRealShift.out GenEigsComplexShift.out \
	SymGEigsCholesky.out SymGEigsRegInv.out

test:
	-./QR.out
	-./Eigen.out
	-./SymEigs.out
	-./SymEigsShift.out
	-./GenEigs.out
	-./GenEigsRealShift.out
	-./GenEigsComplexShift.out
	-./SymGEigsCholesky.out
	-./SymGEigsRegInv.out

%.out: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

clean:
	-rm *.out
