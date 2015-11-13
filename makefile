BUILD := release
cxxflags.debug := -g 
cxxflags.release := -g -O3 -DBDEBUG
CXXFLAGS := ${cxxflags.${BUILD}} -Wall  -std=c++11 -pthread 
INCLUDES := -I /home/sam/programmes/lib/eigen/ -I /home/sam/programmes/lib/boost_1_59_0/
CXX=~/programmes/general/gcc-4.9.2/bin/g++
PROGRAMS := shrek
OBJECTS := usefulTools.o main.o genotypefilehandler.o command.o snp.o region.o genotype.o linkage.o decomposition.o snpestimation.o interval.o 
.PHONY: all clean
shrek: $(OBJECTS)
	$(CXX) $(CXXFLAGS)  $(OBJECTS) -o $@ $(LIBS)
shrek3: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@ $(LIBS)
main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c main.cpp -o main.o
usefulTools.o: usefulTools.cpp
	$(CXX) $(CXXFLAGS) -c usefulTools.cpp -o usefulTools.o
genotypefilehandler.o: genotypefilehandler.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c genotypefilehandler.cpp -o genotypefilehandler.o
command.o: command.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c command.cpp -o command.o
snp.o: snp.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c snp.cpp -o snp.o
region.o: region.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c region.cpp -o region.o
snpestimation.o: snpestimation.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c snpestimation.cpp -o snpestimation.o
genotype.o: genotype.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c genotype.cpp -o genotype.o
linkage.o: linkage.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c linkage.cpp -o linkage.o
interval.o: interval.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c interval.cpp -o interval.o
decomposition.o: decomposition.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c decomposition.cpp -o decomposition.o
all: $(PROGRAMS)
clean:	
	rm -f *.o
