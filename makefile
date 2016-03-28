build := release
cxxflags.debug := -pg -g
cxxflags.release := -O2 -DBDEBUG  -DARMA_NO_DEBUG -mpopcnt
CXXFLAGS := ${cxxflags.${build}} -std=c++11
INCLUDES := -I /home/sam/programmes/lib/armadillo/usr/include/ -I /home/sam/programmes/lib/boost_1_59_0/
LIBS := -larmadillo -pthread -lopenblas -lgomp
CXX=g++
PROGRAMS := shrek
OBJECTS := usefulTools.o main.o command.o snp.o region.o  genotypefilehandler.o snpestimation.o genotype.o linkage.o decomposition.o
.PHONY: all clean
all: $(PROGRAMS)

shrek: $(OBJECTS)
	$(CXX) $(CXXFLAGS)  $(OBJECTS) -o $@ $(LIBS)

main.o: main.cpp command.h region.h region.h genotypefilehandler.h snpestimation.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@ $(LIBS)

usefulTools.o: usefulTools.cpp usefulTools.h
	$(CXX) $(CXXFLAGS) -c  $< -o $@

genotypefilehandler.o: genotypefilehandler.cpp genotypefilehandler.h usefulTools.h snp.h genotype.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c genotypefilehandler.cpp -o genotypefilehandler.o $(LIBS)

command.o: command.cpp command.h usefulTools.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

snp.o: snp.cpp snp.h usefulTools.h region.h command.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

region.o: region.cpp region.h usefulTools.h interval.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

snpestimation.o: snpestimation.cpp snpestimation.h genotypefilehandler.h genotype.h region.h linkage.h decomposition.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

genotype.o: genotype.cpp genotype.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@ $(LIBS)

linkage.o: linkage.cpp linkage.h snp.h genotype.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@ $(LIBS)

decomposition.o: decomposition.cpp decomposition.h linkage.h snp.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@ $(LIBS)
clean:
	rm -f *.o

