/**
 * hamiltonicityMethods.c
 * 
 * Author: Jarne Renders (jarne.renders@kuleuven.be)
 * 
 * A description of the methods can be found in the header file.
 *
 */

#include <stdio.h>
#include <stdbool.h>
#include "bitset.h"
#include "hamiltonicityMethods.h"

bool canBeHamiltonian(bitset adjacencyList[], bitset remainingVertices, int
lastElemOfPath, int firstElemOfPath, int numberOfVertices, int pathLength) {

    // Check whether we have a Hamiltonian path already and whether this path
    // is a cycle.
    if((pathLength == numberOfVertices) && contains(adjacencyList[firstElemOfPath], lastElemOfPath)) {
        return true;
    }

    // Check if cycle can still be closed with remaining vertices.
    if(isEmpty(intersection(adjacencyList[firstElemOfPath],remainingVertices))) { 
        return false;
    }

    // Check for all elements not yet visited whether they still have two
    // neighbours to which they can connect.
    bitset remainingWithFirstAndLast = 
     union(remainingVertices, union(singleton(firstElemOfPath), singleton(lastElemOfPath)));
    forEach(vertex, remainingVertices) {

        //  Neighbours which either do not lie in the path, or which are one
        //  of its endpoints.
        bitset remainingNeighbours = 
         intersection(adjacencyList[vertex], remainingWithFirstAndLast );

        //  If there is only one such neighbour or less, our path cannot be
        //  extended through this vertex into a hamiltonian cycle.
        if(size(remainingNeighbours) < 2) return false;
    }

    // Create a bitset of the neighbours of the last element in the path which
    // do not belong to the path. The path will be extended via these
    // neighbours.
    bitset neighboursOfLastNotInPath = 
     intersection(adjacencyList[lastElemOfPath], remainingVertices);
    forEach(neighbour, neighboursOfLastNotInPath) {

        //  Save the current last element of path. If an extension cannot be a
        //  hamiltonian cycle, we need to try the other possible extensions
        //  for this path.
        int oldElemOfPath = lastElemOfPath;

        //  Extend the path with neighbour, which is a neighbour oldElemOfPath
        //  that does no belong to the path yet.
        removeElement(remainingVertices, neighbour);
        lastElemOfPath = neighbour; // Neighbour is the new last element.

        //  If this extension can become a hamiltonian cycle, so can the
        //  current path.
        if (canBeHamiltonian(adjacencyList, remainingVertices, lastElemOfPath,
         firstElemOfPath, numberOfVertices, pathLength + 1)) {
            return true;
        }

        //  If we reach this part, the extension could not become a
        //  hamiltonian cycle, hence we need to look again at the other
        //  possible extensions for our old path.
        add(remainingVertices, lastElemOfPath);
        lastElemOfPath = oldElemOfPath;
    }

    //  None of the possible extensions worked, so the path cannot be a
    //  hamiltonian cycle.
    return false;
}


bool canBeHamiltonianPrintCycle(bitset adjacencyList[], bitset
remainingVertices, int pathList[], int lastElemOfPath, int firstElemOfPath,
int numberOfVertices, int pathLength, int* numberOfHamiltonianCycles, bool
allCyclesFlag, bool verboseFlag) {

    // Check whether we have a Hamiltonian path already and whether this path is a cycle.
    if((pathLength == numberOfVertices) && contains(adjacencyList[firstElemOfPath], lastElemOfPath)) {
        if(verboseFlag) {
            fprintf(stderr,"Path: ");
            for(int i = 1; i < numberOfVertices; i++) {
                fprintf(stderr, "%d -> ", pathList[i]);
            }
            fprintf(stderr,"%d\n",pathList[0]);
        }
        (*numberOfHamiltonianCycles)++;
        return true;
    }

    // Check if cycle can still be closed with remaining vertices.
    if(isEmpty(intersection(adjacencyList[firstElemOfPath],remainingVertices))) { 
        return false;
    }

    // Check for all elements not yet visited whether they still have two
    // neighbours to which they can connect.
    bitset remainingWithFirstAndLast = 
     union(remainingVertices, union(singleton(firstElemOfPath), singleton(lastElemOfPath)));
    forEach(vertex, remainingVertices) {

        //  Neighbours which either do not lie in the path, or which are one
        //  of its endpoints.
        bitset remainingNeighbours = 
         intersection(adjacencyList[vertex], remainingWithFirstAndLast );

        //  If there is only one such neighbour or less, our path cannot be
        //  extended through this vertex into a hamiltonian cycle.
        if(size(remainingNeighbours) < 2) return false;
    }

    // Create a bitset of the neighbours of the last element in the path which
    // do not belong to the path. The path will be extended via these
    // neighbours.
    bitset neighboursOfLastNotInPath = 
     intersection(adjacencyList[lastElemOfPath], remainingVertices);
    forEach(neighbour, neighboursOfLastNotInPath) {

        //  Save the current last element of path. If an extension cannot be a
        //  hamiltonian cycle, we need to try the other possible extensions
        //  for this path.
        int oldElemOfPath = lastElemOfPath;

        //  Extend the path with neighbour, which is a neighbour oldElemOfPath
        //  that does no belong to the path yet.
        removeElement(remainingVertices, neighbour);
        lastElemOfPath = neighbour; // Neighbour is the new last element.
        pathList[pathLength] = neighbour;

        //  If this extension can become a hamiltonian cycle, so can the
        //  current path.
        if (canBeHamiltonianPrintCycle(adjacencyList, remainingVertices, 
         pathList, lastElemOfPath, firstElemOfPath, numberOfVertices,
         pathLength + 1, numberOfHamiltonianCycles, allCyclesFlag,
         verboseFlag)) {

            // In the case we want to find all cycles, we should only
            // backtrack once we have exhausted all possibilities.
            if(!allCyclesFlag)
                return true;
        }

        //  If we reach this part, the extension could not become a
        //  hamiltonian cycle, hence we need to look again at the other
        //  possible extensions for our old path.
        add(remainingVertices, lastElemOfPath);
        lastElemOfPath = oldElemOfPath;
    }

    //  None of the possible extensions worked, so the path cannot be a
    //  hamiltonian cycle.
    return (*numberOfHamiltonianCycles);
}

bool isHamiltonian(bitset adjacencyList[], int numberOfVertices, bitset
excludedVertices, bool allCyclesFlag, bool verboseFlag) { 
    int numberOfHamiltonianCycles = 0;

    //  We check whether the subgraph spanned by the included vertices is
    //  hamiltonian.
    bitset includedVertices = complement(excludedVertices, numberOfVertices);

    if(isEmpty(includedVertices)) return false;

    // First included vertex.
    int startingVertex = next(includedVertices,-1);
    int lowestDegree = size(adjacencyList[startingVertex]);

    //  Find an included vertex of lowest degree.
    forEachAfterIndex(includedVertex, includedVertices, startingVertex){
        int degreeOfNeighbour = 
         size(intersection(adjacencyList[includedVertex], includedVertices));
        if(lowestDegree > degreeOfNeighbour) {
            lowestDegree = degreeOfNeighbour;
            startingVertex = includedVertex;
        }
    }

    // Loop over included neighbours of startingVertex and for each such
    // neighbour loop over the included neighbours of startingVertex that are
    // of higher index.
    forEach(secondElemOfPath, intersection(adjacencyList[startingVertex], includedVertices)) {
        forEachAfterIndex(lastElemOfPath, intersection(adjacencyList[startingVertex], includedVertices), secondElemOfPath) {

            // Create path, lastElemOfPath, startingVertex, secondElemOfPath.
            // We have lastElemOfPath > secondElemOfPath, so that we
            // eliminate the checking of paths which are mirrored.
            bitset path = singleton(startingVertex);
            add(path, lastElemOfPath);
            add(path, secondElemOfPath);
            bitset remainingVertices = difference(includedVertices, path);

            if(!allCyclesFlag && !verboseFlag) {

                // Check if this path can be extended to some hamiltonian cycle.
                if (canBeHamiltonian(adjacencyList, remainingVertices,
                 lastElemOfPath, secondElemOfPath, size(includedVertices), 3)) {
                    return true;
                }
                continue;
            }

            // If there is a special flag: use canBeHamiltonianPrintCycle
            // (which is a bit slower than canBeHamiltonian).
            int pathList[numberOfVertices];
            pathList[0] = lastElemOfPath;
            pathList[1] = startingVertex;
            pathList[2] = secondElemOfPath;
            canBeHamiltonianPrintCycle(adjacencyList, remainingVertices,
             pathList, secondElemOfPath, lastElemOfPath, size(includedVertices),
             3, &numberOfHamiltonianCycles, allCyclesFlag, verboseFlag);

            //  Stop after one hamiltonian cycle if -a is not present.
            if(!allCyclesFlag && numberOfHamiltonianCycles) {
                return numberOfHamiltonianCycles;
            }
        }
    }
    if(allCyclesFlag) {
       fprintf(stderr,"There were %d hamiltonian cycles in this (sub)graph.\n\n",
        numberOfHamiltonianCycles);
    }

    //  Will be non-zero if there is a hamiltonian cycle.
    return numberOfHamiltonianCycles;
}

bool hasMinimumDegree(bitset adjacencyList[], int numberOfVertices, int
degree) {

    //Check whether all vertices have high enough degree.
    for (int i = 0; i < numberOfVertices; i++) {
        if (size(adjacencyList[i]) < degree) {
            return false;
        }
    }
    return true;
}

bool isK1Hamiltonian(bitset adjacencyList[], int numberOfVertices, bool
verboseFlag, bool allCyclesFlag, int vertexToCheck) {

    //  Graphs with minimum degree < 3 cannot be K1-hamiltonian.
    if(!hasMinimumDegree(adjacencyList,numberOfVertices,3)) {
        if(verboseFlag) {
            fprintf(stderr, "Graph does not have minimum degree 3.\n");
        }
        return false;
    }

    //  An exceptional vertex is one for which the vertex-deleted subgraph is
    //  non-hamiltonian. 
    bitset exceptionalVertices = EMPTY;

    //  Loop over all vertices and determine whether the vertex-deleted
    //  subgraph is hamiltonian.
    for (int i = 0; i < numberOfVertices; i++) {
        bitset excludedVertices = singleton(i);
        if(!verboseFlag) {
            if(!(isHamiltonian(adjacencyList,numberOfVertices,excludedVertices,
             false, false))) {
                return false;
            }
            continue;
        }

        //  The following gets executed only if -v is present.
        bool verbose = false;
        bool cycles = false;

        //  vertexToCheck is determined by -v#.
        if(vertexToCheck == i) {
            verbose = true;
            cycles  = allCyclesFlag;
            fprintf(stderr, "Looking at G - %d.\n", vertexToCheck);
        }
        if (!(isHamiltonian(adjacencyList, numberOfVertices, excludedVertices,
         cycles, verbose))) {
            add(exceptionalVertices, i);
        }
    }

    //  Print out the exceptional vertices.
    int nOfExceptionalVertices = size(exceptionalVertices);
    if(verboseFlag) {
        if(nOfExceptionalVertices) {
            fprintf(stderr, "There are %d exceptional vertices: {",
             nOfExceptionalVertices);
            forEach(excVertex, exceptionalVertices) {
                fprintf(stderr, "%d, ", excVertex);
            }
            fprintf(stderr, "\b\b}\n");
        }
        else {
            fprintf(stderr, "No exceptional vertices.\n");
        }
    }

    //  Zero if there are exceptional vertices, non-zero otherwise.
    return !nOfExceptionalVertices;
}

bool isK2Hamiltonian(bitset adjacencyList[], int numberOfVertices, bool
verboseFlag, bool allCyclesFlag, int vertexPairToCheck[]) {

    //  Graphs with minimum degree < 3 cannot be K2-hamiltonian.
    if(!hasMinimumDegree(adjacencyList,numberOfVertices,3)) {
        if(verboseFlag) {
            fprintf(stderr, "Graph does not have minimum degree 3.\n");
        }
        return false;
    }

    // Pairs (v,w) for which G - v - w is not hamiltonian.
    bitset exceptionalPairs[numberOfVertices];
    if(verboseFlag) {
        for(int i = 0; i < numberOfVertices; i++) {
            exceptionalPairs[i] = EMPTY;
        }
    }
    bool encounteredNonHamSubgraph = false;

    //  Loop over all edges vw with v < w and check if G - v - w is
    //  hamiltonian.
    for (int i = 0; i < numberOfVertices; i++) {
        bitset excludedVertices = singleton(i);
        forEachAfterIndex(neighbour, adjacencyList[i], i) {
            add(excludedVertices, neighbour);
            if(!verboseFlag) {
                if(!(isHamiltonian(adjacencyList, numberOfVertices,
                 excludedVertices, false, false))){
                    return false;
                }
                removeElement(excludedVertices, neighbour);
                continue;
            }

            //  Gets executed if -v is present.
            bool verbose = false;
            bool cycles = false;

            //  vertexPairToCheck is determined by -v#,#
            if((i == vertexPairToCheck[0] && neighbour == vertexPairToCheck[1]) ||
             (i == vertexPairToCheck[1] && neighbour == vertexPairToCheck[0])) {
                verbose = true;
                cycles = allCyclesFlag;
                fprintf(stderr, "Looking at G - %d - %d.\n",
                 vertexPairToCheck[0], vertexPairToCheck[1]);
            }
            if(!(isHamiltonian(adjacencyList, numberOfVertices,
             excludedVertices, cycles, verbose))){
                add(exceptionalPairs[i], neighbour);
                encounteredNonHamSubgraph = true;
            }
            removeElement(excludedVertices, neighbour);
        }
    }

    if(verboseFlag) {
        if(encounteredNonHamSubgraph) {
            fprintf(stderr, "G - v - w is not hamiltonian for (v,w) in {");
            for (int v = 0; v < numberOfVertices; v++) {
                forEachAfterIndex(w, exceptionalPairs[v], v) {
                    fprintf(stderr, "(%d,%d), ",v,w);
                }
            }
            fprintf(stderr, "\b\b}\n");
        }
        else {
            fprintf(stderr, "Graph is K2-hamiltonian.\n");
        }
    }
    return !encounteredNonHamSubgraph;
}

int containsHamiltonianPathWithEnds(bitset adjacencyList[], int
numberOfVertices, bitset excludedVertices, int start, int end, bool
allCyclesFlag, bool verboseFlag) {
    
    //  If start or end are excluded there cannot be a path between them.
    if(contains(excludedVertices,start) || contains(excludedVertices,end)) {
        return false;
    }

    bitset path = union(singleton(start), singleton(end));
    bitset includedVertices = complement(excludedVertices, numberOfVertices);
    bitset remainingVertices = difference(includedVertices, path);
    if(!verboseFlag && !allCyclesFlag) {

        //  Will return true if this path can be extended to a hamiltonian
        //  path between start and end and false otherwise..
        return canBeHamiltonian(adjacencyList, remainingVertices, start, end,
         size(includedVertices), 2);
    }

    //  Only gets executed if -v or -a are present.
    int pathList[size(includedVertices)];
    pathList[0] = end;
    pathList[1] = start;
    int nOfPaths = 0;

    canBeHamiltonianPrintCycle(adjacencyList, remainingVertices, pathList,
    start, end, size(includedVertices), 2, &nOfPaths, allCyclesFlag,
    verboseFlag);

    if(allCyclesFlag) {
       fprintf(stderr,"There were %d hamiltonian (%d,%d)-paths in this graph.\n\n",
        nOfPaths, start, end);
    }

    //  Will return 0 if there are no hamiltonian paths and non-zero
    //  if there are.
    return nOfPaths;
}

bool isPartOfDisjointSpanningPaths(bitset adjacencyList[], bitset currentPath,
bitset excludedVertices, int pathList[], int firstElemOfPath, int
lastElemOfPath, bitset verticesContainedByPath1, int numberOfVertices, int
firstElemOfPath2, int lastElemOfPath2, bitset verticesContainedByPath2, int*
nOfSpanningPaths, bool allCyclesFlag, bool verboseFlag)  { 

    //  Check for second path if first path is a cycle and contains
    //  all required vertices.
    if(contains(adjacencyList[lastElemOfPath], firstElemOfPath) &&
     equals(intersection(currentPath, verticesContainedByPath1), verticesContainedByPath1)) {

        //  Included vertices which do not belong to the first path
        //  and which are not the endpoints of the second path.
        bitset remainingVertices = complement(union(currentPath,excludedVertices),
         numberOfVertices+size(excludedVertices));
        removeElement(remainingVertices, firstElemOfPath2);
        removeElement(remainingVertices, lastElemOfPath2); 
        if(verboseFlag || allCyclesFlag) {
            int secondPath[numberOfVertices];
            secondPath[0] = firstElemOfPath2;
            secondPath[1] = lastElemOfPath2;
            int nOfPaths = 0;

            //  Check whether the subgraph spanned by the remaining vertices
            //  contains a hamiltonian path between firstElemOfPath2 and
            //  lastElemOfPath2.
            if(canBeHamiltonianPrintCycle(adjacencyList, remainingVertices,
             secondPath, lastElemOfPath2, firstElemOfPath2,
             size(remainingVertices)+2, 2, &nOfPaths, allCyclesFlag,
             verboseFlag)) {
                (*nOfSpanningPaths) += nOfPaths;
                if(verboseFlag) {fprintf(stderr,"Second path: ");
                    for(int i = 1; i < size(currentPath); i++) {
                        fprintf(stderr, "%d -> ", pathList[i]);
                    }
                    fprintf(stderr,"%d\n",pathList[0]);
                }
                if(!allCyclesFlag) {
                    return true;
                }
            }
        }
        else {

            //  Check whether the subgraph spanned by the remaining vertices
            //  contains a hamiltonian path between firstElemOfPath2 and
            //  lastElemOfPath2.
            if(canBeHamiltonian(adjacencyList, remainingVertices,
             lastElemOfPath2, firstElemOfPath2,
             size(remainingVertices)+2, 2)) {
                return true;
            }
        }
    }

    //  Included neighbours of the last element of path1 which are
    //  themselves not in path1 and which are not equal to the ends
    //  of path2.
    bitset remainingNeighboursOfLast = difference(adjacencyList[lastElemOfPath],
     union(excludedVertices,union(currentPath, union(verticesContainedByPath2,
     union(singleton(firstElemOfPath2), singleton(lastElemOfPath2))))));
    forEach(neighbour, remainingNeighboursOfLast) {
        pathList[size(currentPath)] = neighbour;
        add(currentPath, neighbour);
        if(isPartOfDisjointSpanningPaths(adjacencyList, currentPath, 
         excludedVertices, pathList, firstElemOfPath, neighbour, 
         verticesContainedByPath1, numberOfVertices, firstElemOfPath2,
         lastElemOfPath2, verticesContainedByPath2, nOfSpanningPaths,
         allCyclesFlag, verboseFlag)) {
            if(!allCyclesFlag) {
                return true;
            }
        }
        removeElement(currentPath, neighbour);
    }
    return (*nOfSpanningPaths);
};


bool containsDisjointSpanningPathsWithEnds(bitset adjacencyList[], int
numberOfVertices, bitset excludedVertices, int startOfPath1, int endOfPath1,
bitset verticesContainedByPath1, int startOfPath2, int endOfPath2, bitset
verticesContainedByPath2, bool allCyclesFlag, bool verboseFlag) {

    int nOfPaths = 0;
    bitset path1 = union(singleton(startOfPath1), singleton(endOfPath1));
    int path1List[numberOfVertices - size(excludedVertices)];
    path1List[0] = endOfPath1;
    path1List[1] = startOfPath1;

    bool isPart = isPartOfDisjointSpanningPaths(adjacencyList, path1,
    excludedVertices, path1List, endOfPath1, startOfPath1,
    verticesContainedByPath1, numberOfVertices - size
    (excludedVertices), startOfPath2, endOfPath2,
    verticesContainedByPath2, &nOfPaths, allCyclesFlag,
    verboseFlag); 

    if(allCyclesFlag) {
        fprintf(stderr, "Graph contains %d pairs of disjoint spanning paths between (%d,%d)%s and (%d,%d)%s.\n",
         nOfPaths,startOfPath1, endOfPath1,
         isEmpty(verticesContainedByPath1) ? "" : " containing specified vertices",
         startOfPath2, endOfPath2,
         isEmpty(verticesContainedByPath1) ? "" : " containing specified vertices");
    }

    return isPart;
}

bool isTraceable(bitset adjacencyList[], int numberOfVertices, bitset
excludedVertices, bool allCyclesFlag, bool verboseFlag) {
    long long unsigned nOfPaths = 0;
    for(int i = 0; i < numberOfVertices; i++) {
        for(int j = i + 1; j < numberOfVertices; j++) {
            int nOfPathsWithEnds;
            if((nOfPathsWithEnds = containsHamiltonianPathWithEnds(
             adjacencyList, numberOfVertices, excludedVertices, i, j,
             allCyclesFlag, verboseFlag))) {
                if(!allCyclesFlag) {
                    return true;
                }
                nOfPaths += nOfPathsWithEnds;
            }
        }
    }
    if(allCyclesFlag) {
        fprintf(stderr, "There were %llu hamiltonian paths in the (sub)graph\n", nOfPaths);
    }
    return nOfPaths;
}

bool isK1Traceable(bitset adjacencyList[], int numberOfVertices, bool
allCyclesFlag, bool verboseFlag, int vertexToCheck) {

    //  An exceptional vertex is one for which the vertex-deleted subgraph is
    //  non-traceable. 
    bitset exceptionalVertices = EMPTY;

    //  Loop over all vertices and determine whether the vertex-deleted
    //  subgraph is traceable.
    for (int i = 0; i < numberOfVertices; i++) {
        bitset excludedVertices = singleton(i);
        if(!verboseFlag) {
            if(!(isTraceable(adjacencyList,numberOfVertices,excludedVertices,
             false, false))) {
                return false;
            }
            continue;
        }

        //  The following gets executed only if -v is present.
        bool verbose = false;
        bool cycles = false;

        //  vertexToCheck is determined by -v#.
        if(vertexToCheck == i) {
            verbose = true;
            cycles  = allCyclesFlag;
            fprintf(stderr, "Looking at G - %d.\n", vertexToCheck);
        }
        if (!(isTraceable(adjacencyList, numberOfVertices, excludedVertices,
         cycles, verbose))) {
            add(exceptionalVertices, i);
        }
    }

    //  Print out the exceptional vertices.
    int nOfExceptionalVertices = size(exceptionalVertices);
    if(verboseFlag) {
        if(nOfExceptionalVertices) {
            fprintf(stderr, "There are %d exceptional vertices: {",
             nOfExceptionalVertices);
            forEach(excVertex, exceptionalVertices) {
                fprintf(stderr, "%d, ", excVertex);
            }
            fprintf(stderr, "\b\b}\n");
        }
        else {
            fprintf(stderr, "No exceptional vertices.\n");
        }
    }

    //  Zero if there are exceptional vertices, non-zero otherwise.
    return !nOfExceptionalVertices;
}
