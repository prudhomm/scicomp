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

SOURCES=test_jacobi.cpp test_legendreset.cpp
test_jacobi_OBJS=test_jacobi.o
test_legendreset_OBJS=test_legendreset.o

all: test_jacobi test_legendreset


test_jacobi: $(test_jacobi_OBJS)
	$(CXX) $(LDFLAGS) -o test_jacobi $(test_jacobi_OBJS) $(LIBS)

test_legendreset: $(test_legendreset_OBJS)
	$(CXX) $(LDFLAGS) -o test_legendreset $(test_legendreset_OBJS) $(LIBS)

.cpp.o: $(SOURCES)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

.SUFFIXES: .cpp .o

test_jacobi.o: ublas.hpp jacobi.hpp
test_legendreset.o: ublas.hpp jacobi.hpp legendreset.hpp

### Makefile ends here
