all: build


GITVERSION= "no-version"

CXXFLAGS+= -std=gnu++11 -DGITVERSION=\"$(GITVERSION)\"
CXXFLAGS+= -Wno-format -Wno-pointer-arith -I ../

CXXFLAGS+=   -march=native -DHASHSEED1=0x74b0dc51 -DHASHSEED2=0x19495cff -g -D VALUELENGTH=12


HEADERS= ../*.h

build: *.cc


run: build
	./build overall.txt /home/yy/gitothello/Othello/src/shortbin/ .Kmer.bin.short 
