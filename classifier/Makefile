all:

OPTFLAGS = -O0 -g -rdynamic -std=c++0x
OUTPUTDIR = ./

CFLAGS += $(OPTFLAGS) -fopenmp

#CXXFLAGS= -std=c++11 -D USEBYTE -D NEWHASH -D PREFETCH -g 

GITVERSION= "no-version"

CXXFLAGS+= -std=gnu++11 -DGITVERSION=\"$(GITVERSION)\"
CXXFLAGS+= -Wno-format -Wno-pointer-arith
CXXFLAGS+= -DKMERLENGTH=20   -march=native -DHASHSEED1=0x74b0dc51 -DHASHSEED2=0x19495cff -O3

CXXFLAGS+= -fopenmp

#SRC= src/othello/testOthelloKeyValue.cc

HEADERS= src/othello/*.h

classifier: src/metagenomics/assignMetagenomicsRead_allTaxoRank.cpp src/othello/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignMetagenomicsRead_allTaxoRank.cpp -DVALUELENGTH=12 -D W2

assignMetagenomicsRead_allTaxoRank_12_w2: src/metagenomics/assignMetagenomicsRead_allTaxoRank.cpp src/othello/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignMetagenomicsRead_allTaxoRank.cpp -DVALUELENGTH=12 -D W2

assignMetagenomicsRead_allTaxoRank_13_w2: src/metagenomics/assignMetagenomicsRead_allTaxoRank.cpp src/othello/*.h
	g++ $(CXXFLAGS) -o $@ src/metagenomics/assignMetagenomicsRead_allTaxoRank.cpp -DVALUELENGTH=13 -D W2

getTaxonReadCount: 
	g++ $(CXXFLAGS) -o $@ src/others/utils/getTaxonReadCount.cpp
	