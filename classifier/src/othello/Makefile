all:	SortAndCompressKV	testMulIndex testRatio test4


#GITVERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
GITVERSION= "no-version"

CXXFLAGS+= -std=gnu++11 -DGITVERSION=\"$(GITVERSION)\"
CXXFLAGS+= -Wno-format -Wno-pointer-arith

CXXFLAGS+= -DKMERLENGTH=20   -march=native -DHASHSEED1=0x74b0dc51 -DHASHSEED2=0x19495cff -g

#CXXFLAGS+= -DKMERLENGTH=20 -march=native -g


SRC= testOthelloKeyValue.cc

HEADERS= *.h

.PHONY: clean

test2: $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC) -DVALUELENGTH=2

test4: $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC) -DVALUELENGTH=4

test8: $(SRC) $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $(SRC) -DVALUELENGTH=8

test6: $(SRC) $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $(SRC) -DVALUELENGTH=6

testLargeOthello: testLargeOthello.cc $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< -DVALUELENGTH=6

test5: $(SRC) $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $(SRC) -DVALUELENGTH=5

test1: $(SRC) $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $(SRC) -DVALUELENGTH=1

test16: $(SRC) $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $(SRC) -DVALUELENGTH=16

testMulIndex: testMulIndex.cc $(HEADERS)
	g++ $(CXXFLAGS) -o $@ $< -DVALUELENGTH=16

testIndex: testIndex.cc *.h
	g++ $(CXXFLAGS) -o $@ testIndex.cc

clean:	
	rm test

run: test4 test1 test2 test5
	./test5 1 smallkmers.txt tmp
	./test4 2 smallkmers.txt tmp
	./test2 2 smallkmers.txt tmp
	./test1 1 smallkmers.txt tmp

testRatio: testRatio.cc othello.h

doc:
	doxygen Doxyfile
