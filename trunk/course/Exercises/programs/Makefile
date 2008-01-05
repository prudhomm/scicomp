CXX=g++
CXXFLAGS=-O2
CPPFLAGS=-I$(HOME)/blitz/include -I$(HOME)/gmm/include
LDFLAGS=-L$(HOME)/blitz/lib 
LIBS=-lblitz -lblas
PROGRAMS=td1_exo1 td1_exo2
all: $(PROGRAMS)
# link
$(PROGRAMS): 
	$(CXX) -o $@ $@.o $(LDFLAGS) $(LIBS)
# compile
.cpp.o: 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS)  -c $*.cpp
# clean the directory
clean: 
	-rm -f *.o $(TARGETS)
# dependencies
td1_exo1:      td1_exo1.o 
td1_exo2:      td1_exo2.o 

.SUFFIXES: .cpp .o .hpp
