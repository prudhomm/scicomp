### Makefile ---

## Author: christophe.prudhomme@ujf-grenoble.fr
## Version: $Id: Makefile,v 0.0 2008/01/25 17:20:48 prudhomm Exp $
## Keywords:
## X-URL:

CXX=g++
CXXFLAGS=-O2
LDFLAGS=
LIBS=-lm

CPPFLAGS=-I.

SOURCES=test_jacobi.cpp
OBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: test_jacobi


test_jacobi: $(OBJS)
	$(CXX) $(LDFLAGS) -o test_jacobi $(OBJS) $(LIBS)

.cpp.o: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

.SUFFIXES: .cpp .o

test_jacobi.o: ublas.hpp jacobi.hpp

### Makefile ends here
