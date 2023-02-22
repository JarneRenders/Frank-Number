# FindFrankNumber
This repository contains a program created for the article "J. Goedgebeur, E. Máčajová, J. Renders: Frank number and nowhere-zero flows on
graphs, manuscript".

The program uses Brendan McKay's graph6 format to read and write graphs. See <http://users.cecs.anu.edu.au/~bdm/data/formats.txt>.

### Short manual
This program can be used to determine whether a given 3-edge-connected cubic graph has Frank number 2 or not, however, without any optional parameters it is assumed all input graphs are cyclically 4-edge-connected cubic graphs. The program makes use of two algorithms, a heuristic algorithm which checks sufficient conditions for graphs to have Frank number 2 and an exact algorithm. This heuristic algorithm only works for cyclically 4-edge-connected graphs. Without any extra flags first the sufficient condition is test and if it fails the exact algorithm is performed. 

This program supports graphs up to 128 vertices.

### Installation

This requires a working shell and `make`. Navigate to the folder containing findFrankNumber.c and compile using:

* `make` to create a binary for the 64-bit version;
* `make 128bit` to create a binary for the 128-bit version;
* `make 128bitarray` to create a binary for an alternative 128-bit version;
* `make all` to create all the above binaries.

The 64-bit version supports graphs only up to 64 vertices, the 128-bit versions up to 128 vertices. For graphs containing up to 64 vertices the 64-bit version performs siginificantly faster than the 128-bit versions. Typically, the 128-bit array version performs faster than the standard 128-bit version. Use `make clean` to remove all binaries created in this way.

### Usage of findFrankNumber

All options can be found by executing `./findFrankNumber -h`.

Usage: `./findFrankNumber [-2|-e] [-b] [-c] [-d] [-h] [-p] [-s] [-v] [res/mod]`

Filter 3-edge-connected cubic graphs having Frank number 2. Unless option -e is present, correct output is only guaranteed if the graphs are also cyclically 4-edge-connected. By default, an input graph will be send to stdout if its Frank number is not equal to 2.

Graphs are read from stdin in graph6 format. Graphs are sent to stdout in graph6 format. If the input graph had a graph6 header, so will the output graph (if it passes through the filter).

The order in which the arguments appear does not matter.
```
  -2, --only-heuristic          Only perform the heuristic algorithm, i.e.
                                 check whether the graph passes the sufficient
                                 condition; The heuristic algorithm only works
                                 for cyclically 4-edge-connected graphs
  -b, --brute-force             Whenever a graph is checked using the exact 
                                 algorithm apply a brute force method instead
  -c, --complement              Reverse output of the graphs, i.e. output all 
                                 graphs which would not be output without this
                                 flag and do not output those which would
  -d, --double-check            Whenever a graph passes the sufficient
                                 condition, double check the result by 
                                 computing the corresponding orientations
  -e, --only-exact              Only perform the exact algorithm and not the 
                                 heuristic one; This flag needs to be present
                                 for graphs which are not cyclically 
                                 4-edge-connected
  -h, --help                    Print this help text
  -p, --print-orientation       Print the two orientations for graphs 
                                 determined to have Frank number 2
  -s, --single-graph-parallel   Parallellize the computation of the exact
                                 method for a single graph; Use with res/mod
  -v, --verbose                 Give more detailed output
  res/mod                       Split the generation in mod (not necessarily
                                 equally big) parts; Here part res will be 
                                 executed
```

### Examples

`./findFrankNumber`
Sends all graphs for which the Frank number is not equal to 2 from stdin to stdout. Correct output is only guaranteed for cyclically 4-edge-connected cubic graphs.

`./findFrankNumber -c`
Sends all graphs for which the Frank number is equal to 2 from stdin to stdout. Correct output is only guaranteed for cyclically 4-edge-connected cubic graphs.

`./findFrankNumber -2`
Sends all graphs for which the sufficient condition fails from stdin to stdout. Correct output is only guaranteed for cyclically 4-edge-connected cubic graphs.

`./findFrankNumber -e`
Sends all graphs for which the Frank number is not equal to 2 from stdin to stdout. Correct output is only guaranteed for 3-edge-connected cubic graphs.

`./findFrankNumber -p`
Prints to stderr for all graphs which are determined to have Frank number 2, the two orientations showing this. Correct output is only guaranteed for cyclically 4-edge-connected cubic graphs.

`./findFrankNumber 3/8`
The same behaviour as `./findFrankNumber`, but only consider every eigth graph from stdin starting from the third graph.

`./findFrankNumber -s 3/8`
The same behaviour as `./findFrankNumber`, but the computation is parallellized for a single graph. This does not parallelize the heuristic algorithm, but it does for the exact algorithm. If some part determines that the Frank number is 2, this is the case. Only if all parts cannot determine the Frank number is 2, the Frank number is not equal to 2.