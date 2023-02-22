compiler=gcc
flags=-std=gnu11 -march=native -Wall -Wno-missing-braces

# The 64-bit version of this program is faster but only supports graphs up to 64 vertices.
64bit: findFrankNumber.c readGraph/readGraph6.c bitset.h hamiltonicityMethods.c
	$(compiler) -DUSE_64_BIT -o findFrankNumber findFrankNumber.c readGraph/readGraph6.c hamiltonicityMethods.c $(flags) -O3

128bit: findFrankNumber.c readGraph/readGraph6.c bitset.h hamiltonicityMethods.c
	$(compiler) -DUSE_128_BIT -o findFrankNumber-128 findFrankNumber.c readGraph/readGraph6.c hamiltonicityMethods.c $(flags) -O3

128bitarray: findFrankNumber.c readGraph/readGraph6.c bitset.h hamiltonicityMethods.c
	$(compiler) -DUSE_128_BIT_ARRAY -o findFrankNumber-128a findFrankNumber.c readGraph/readGraph6.c hamiltonicityMethods.c $(flags) -O3

profile: findFrankNumber.c readGraph/readGraph6.c bitset.h hamiltonicityMethods.c
	$(compiler) -DUSE_64_BIT -o findFrankNumber-pr findFrankNumber.c readGraph/readGraph6.c hamiltonicityMethods.c $(flags) $(densenauty32) -g -pg

all: 64bit 128bit 128bitarray

.PHONY: clean
clean:
	rm -f findFrankNumber findFrankNumber-128 findFrankNumber-128a findFrankNumber-pr

