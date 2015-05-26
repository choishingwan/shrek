BUILD := debug
cxxflags.debug := -g 
cxxflags.release := -g -O3 -DBDEBUG
CXXFLAGS := ${cxxflags.${BUILD}} -Wall  -std=c++0x -pthread
INCLUDES := -I /home/sam/programmes/lib/eigen/
CXX=g++
PROGRAMS := shrek
OBJECTS := usefulTools.o snpindex.o main.o genotypefilehandler.o command.o snp.o region.o decomposition.o decompositionthread.o genotype.o linkage.o linkagethread.o snpestimation.o
.PHONY: all clean
shrek: $(OBJECTS)
	$(CXX) $(CXXFLAGS)  $(OBJECTS) -o $@ $(LIBS)
main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c main.cpp -o main.o
usefulTools.o: usefulTools.cpp
	$(CXX) $(CXXFLAGS) -c usefulTools.cpp -o usefulTools.o
snpindex.o: snpindex.cpp
	$(CXX) $(CXXFLAGS) -c snpindex.cpp -o snpindex.o 
genotypefilehandler.o: genotypefilehandler.cpp
	$(CXX) $(CXXFLAGS) -c genotypefilehandler.cpp -o genotypefilehandler.o
command.o: command.cpp
	$(CXX) $(CXXFLAGS) -c command.cpp -o command.o
snp.o: snp.cpp
	$(CXX) $(CXXFLAGS) -c snp.cpp -o snp.o
region.o: region.cpp
	$(CXX) $(CXXFLAGS) -c region.cpp -o region.o
snpestimation.o: snpestimation.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c snpestimation.cpp -o snpestimation.o
decomposition.o: decomposition.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c decomposition.cpp -o decomposition.o
decompositionthread.o: decompositionthread.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c decompositionthread.cpp -o decompositionthread.o
genotype.o: genotype.cpp
	$(CXX) $(CXXFLAGS) -c genotype.cpp -o genotype.o
linkage.o: linkage.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c linkage.cpp -o linkage.o
linkagethread.o: linkagethread.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c linkagethread.cpp -o linkagethread.o
all: $(PROGRAMS)
clean:	
	rm -f $(PROGRAMS) *.o
