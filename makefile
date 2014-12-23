BUILD := debug
cxxflags.debug := -g 
cxxflags.release := -g -O3 -DBDEBUG
CXXFLAGS := ${cxxflags.${BUILD}} -Wall  -std=c++0x
INCLUDES := -I /home/sam/programmes/lib/eigen/
LIBS :=
CXX=g++
PROGRAMS := shrek
OBJECTS := usefulTools.o snpindex.o main.o genotypefilehandler.o command.o snp.o region.o snpestimation.o
.PHONY: all clean
shrek: $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LIBS)
%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(input) -o $(output)

all: $(PROGRAMS)
clean:	
	rm -f $(PROGRAMS) *.o
