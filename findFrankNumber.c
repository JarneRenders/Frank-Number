/**
 * findFrankNumber.c
 * 
 * Author: Jarne Renders (jarne.renders@kuleuven.be)
 *
 */

#define USAGE \
"\nUsage: `./findFrankNumber [-2|-e] [-b] [-c] [-d] [-h] [-p] [-s] [-v] [res/mod]`\n"
#define HELPTEXT \
"Filter 3-edge-connected cubic graphs having Frank number 2.\n\
Unless option -e is present, correct output is only guaranteed if the graphs\n\
are also cyclically 4-edge-connected. By default, an input graph will be send\n\
to stdout if its Frank number is not equal to 2.\n\
\n\
Graphs are read from stdin in graph6 format. Graphs are sent to stdout in\n\
graph6 format. If the input graph had a graph6 header, so will the output\n\
graph (if it passes through the filter).\n\
\n\
The order in which the arguments appear does not matter.\n\
\n\
  -2, --only-heuristic          Only perform the heuristic algorithm, i.e.\n\
                                 check whether the graph passes the sufficient\n\
                                 condition; The heuristic algorithm only works\n\
                                 for cyclically 4-edge-connected graphs\n\
  -b, --brute-force             Whenever a graph is checked using the exact\n\
                                 algorithm apply a brute force method instead\n\
  -c, --complement              Reverse output of the graphs, i.e. output all\n\
                                 graphs which would not be output without this\n\
                                 flag and do not output those which would\n\
  -d, --double-check            Whenever a graph passes the sufficient\n\
                                 condition, double check the result by\n\
                                 computing the corresponding orientations\n\
  -e, --only-exact              Only perform the exact algorithm and not the\n\
                                 heuristic one; This flag needs to be present\n\
                                 for graphs which are not cyclically\n\
                                 4-edge-connected\n\
  -h, --help                    Print this help text\n\
  -p, --print-orientation       Print the two orientations for graphs\n\
                                 determined to have Frank number 2\n\
  -s, --single-graph-parallel   Parallellize the computation of the exact\n\
                                 method for a single graph; Use with res/mod\n\
  -v, --verbose                 Give more detailed output\n\
  res/mod                       Split the generation in mod (not necessarily\n\
                                 equally big) parts; Here part res will be\n\
                                 executed\n\
"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <string.h>
#include "readGraph/readGraph6.h"
#include "bitset.h"
#include "hamiltonicityMethods.h"

struct counters {
    long long unsigned int generatedOrientations;
    long long unsigned int mostGeneratedOrientations;
    long long unsigned int storedBitsets;
    long long unsigned int mostStoredBitsets;
    long long unsigned int orientationsGivingSubset;
    long long unsigned int orientationsGivingSuperset;
    long long unsigned int emptyBitsetsStored;
    long long unsigned int complementaryBitsets;
    long long unsigned int graphsSatisfyingOddnessCondition;
    long long unsigned int graphsNotSatisfyingOddnessCondition;
    long long unsigned int graphsSatisfyingFirstOddness;
    long long unsigned int graphsSatisfyingSecondOddness;
    long long unsigned int totalOrientationsGenerated;
};

struct options {
    bool bruteForceFlag;
    bool complementFlag;
    bool doublecheckFlag;
    bool exhaustiveCheckFlag;
    bool oddCyclesHeuristicFlag;
    bool verboseFlag;
    bool printFlag;
    bool singleGraphFlag;
    int modulo;
    int remainder;
    unsigned long long int sizeOfArray;
};

//******************************************************************************
//
//                          Dynamic arrays
//
//******************************************************************************

typedef struct {
  bitset *array;
  size_t used;
  size_t size;
} Array;

void initArray(Array *a, size_t initialSize) {
  a->array = malloc(initialSize * sizeof(bitset));
  if(a->array == NULL) {
    fprintf(stderr, "Error: out of memory\n");
    exit(1);
  }
  a->used = 0;
  a->size = initialSize;
}

// Double arraysize when too big.
void insertArray(Array *a, bitset element) {
  if (a->used == a->size) {
    a->size *= 2;
    a->array = realloc(a->array, a->size * sizeof(bitset));
    if(a->array == NULL) {
        fprintf(stderr, "Error: out of memory\n");
        exit(1);
    }
  }
  a->array[a->used++] = element;
}

void insertArrayAtPos(Array *a, bitset element, size_t index) {
    if (index > a->size) {
        fprintf(stderr, "Error: index does not lie in the array.\n");
        exit(1);
    }
  a->array[index] = element;
}

void freeArray(Array *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}

//******************************************************************************
//
//                          Digraphs
//
//******************************************************************************

// Digraph structure containing bitset representations.
struct diGraph {
    int numberOfVertices;
    bitset* adjacencyList;
    bitset* reverseAdjacencyList;
    int numberOfArcs;
}; 

//  Initializer for empty graph.
#define emptyGraph(g) {\
 (g)->numberOfArcs = 0;\
 for(int i = 0; i < (g)->numberOfVertices; i++) {\
    (g)->adjacencyList[i] = EMPTY;\
    (g)->reverseAdjacencyList[i] = EMPTY;\
 }\
}

//  Add one directed edge. numberOfArcs will be incorrect if adding existing
//  arc.
#define addArc(g,i,j) {\
 add((g)->adjacencyList[i], j); (g)->numberOfArcs++;\
 add((g)->reverseAdjacencyList[j], i);\
}

//  Remove one undirected edge. numberOfArcs will be incorrect if removing
//  non-existing arc.
#define removeArc(g,i,j) {\
 removeElement((g)->adjacencyList[i], j);(g)->numberOfArcs--;\
 removeElement((g)->reverseAdjacencyList[j], i);\
}

//  Print adjacency list of digraph.
void printDiGraph(struct diGraph *g) {
    for(int i = 0; i < g->numberOfVertices; i++) {
        fprintf(stderr, "%d:", i);
        forEach(nbr, g->adjacencyList[i]) {
            fprintf(stderr, " %d", nbr);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr,"\n");   
}

//  Print adjacency list.
void printGraph(bitset adjacencyList[], int numberOfVertices) {
    for(int i = 0; i < numberOfVertices; i++) {
        fprintf(stderr, "%d: ", i);
        forEach(nbr, adjacencyList[i]) {
            fprintf(stderr, "%d ", nbr);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

//******************************************************************************
//
//                         Strong connectivity check
//
//******************************************************************************

void visit(struct diGraph *g, int vertexToVisit, bitset *unvisitedVertices,
 int L[], int *lengthOfL) {
    if(!contains(*unvisitedVertices, vertexToVisit)) {
        return;
    }
    removeElement(*unvisitedVertices, vertexToVisit);
    forEach(outNeighbour, g->adjacencyList[vertexToVisit]) {
        visit(g, outNeighbour, unvisitedVertices, L, lengthOfL);
    }
    L[*lengthOfL] = vertexToVisit;
    (*lengthOfL)++;
}

void assign(struct diGraph *g, int vertex, bitset *assignedVertices) {
    if(contains(*assignedVertices, vertex)) {
        return;
    }
    add(*assignedVertices, vertex);
    forEach(inNeighbour, g->reverseAdjacencyList[vertex]) {
        assign(g, inNeighbour, assignedVertices);
    }
}

bool isStronglyConnected(struct diGraph *g) {
    bitset unvisitedVertices = complement(EMPTY, g->numberOfVertices);
    int L[g->numberOfVertices];
    int lengthOfL = 0;
    for(int i = 0; i < g->numberOfVertices ; i++) {
        visit(g, i, &unvisitedVertices, L, &lengthOfL);
    }
    bitset assignedVertices = EMPTY;
    assign(g, L[lengthOfL-1], &assignedVertices);

    return size(assignedVertices) == g->numberOfVertices;
}
//******************************************************************************
//
//                     Deletable edges
//
//******************************************************************************

//  Give all edges a graph an index from 0 to |E(G)| - 1 and store in matrix
//  edgeIndices.
void numberEdges(bitset adjacencyList[], int numberOfVertices,
 int edgeIndices[][numberOfVertices]) {
    int counter = 0;
    for(int i = 0; i < numberOfVertices; i++) {
        forEachAfterIndex(nbr, adjacencyList[i], i) {
            edgeIndices[i][nbr] = counter;
            edgeIndices[nbr][i] = counter;
            counter++;
        }
    }
}

//  Used for checking if edge is deletable.
bool containsDirectedPathBetween(struct diGraph *orientation,
 bitset unvisitedVertices, int i, int end) {

    if(contains(orientation->adjacencyList[i], end)) {
        return true;
    }
    removeElement(unvisitedVertices, i);
    forEach(element, intersection(orientation->adjacencyList[i],
     unvisitedVertices)) {
        if(containsDirectedPathBetween(orientation, unvisitedVertices, element,
         end)) {
            return true;
        }
    }
    return false;
}

//  We assume that the given orientation is strongly connected.
bitset getDeletableEdges(struct diGraph *orientation, int numberOfVertices,
 int edgeNumbering[][numberOfVertices]) {

    bitset deletableEdges = EMPTY;

    for(int i = 0; i < numberOfVertices; i++) {
        forEach(nbr, orientation->adjacencyList[i]) {
            removeArc(orientation, i, nbr);
            if(containsDirectedPathBetween(orientation,
             complement(EMPTY, numberOfVertices), i, nbr)) {
                add(deletableEdges, edgeNumbering[i][nbr]);
            }
            addArc(orientation, i, nbr);
        }
    }

    return deletableEdges;
}

void printDeletableEdges(int numberOfVertices,
 int edgeNumbering[][numberOfVertices], bitset orientation[], 
 bitset deletableEdges) {
    fprintf(stderr, "Deletable edges: ");
    for(int i = 0; i < numberOfVertices; i++) {
        forEach(nbr, orientation[i]) {
            if(contains(deletableEdges,edgeNumbering[i][nbr])) {
                fprintf(stderr, "%d--%d ", i, nbr);
            }
        }
    }
    fprintf(stderr, "\n");
}

//******************************************************************************
//
//                          Exact algorithm
//
//******************************************************************************

#define isSubset(set1, set2) equals((set1), intersection((set1),(set2))) 

// Brute force approach
int getIntermediateFrankNumber(struct options *options,
 struct counters *numberOf, int numberOfVertices,
 int edgeNumbering[][numberOfVertices], Array *bitsetsOfDeletableEdges,
 bitset deletableEdges) {

    size_t insertPosition = bitsetsOfDeletableEdges->used;
    bitset bitsetContainingAllEdges = complement(EMPTY, 3*numberOfVertices/2);
    bitset *array = bitsetsOfDeletableEdges->array;

    // Check if Frank number is 2
    for(size_t i = 0; i < bitsetsOfDeletableEdges->used; i++) {

        if(!isEmpty(array[i])) {

            //  If the deletable edges of new orientation is a subset of older
            //  we can dismiss it.
            if(isSubset(deletableEdges, array[i])) {
                numberOf->orientationsGivingSubset++;
                return 0;
            }

            //  If the deletable edges of new orientation is superset of older
            //  we can dismiss older. We set it to EMPTY.
            if(isSubset(array[i], deletableEdges)) {
                if(insertPosition == bitsetsOfDeletableEdges->used) {
                    numberOf->orientationsGivingSuperset++;
                }
                array[i] = EMPTY;
            }

            //  If union of new and older deletable edges are all edges, Frank
            //  number is 2.
            if(equals(union(deletableEdges, array[i]),
             bitsetContainingAllEdges)) {
                numberOf->complementaryBitsets++;
                insertArray(bitsetsOfDeletableEdges, deletableEdges);
                return 2;
            }
        }
        else {

            //  Prepare to store new deletable edges in first EMPTY position.
            if(insertPosition == bitsetsOfDeletableEdges->used) {
                insertPosition = i;
            }
        }
    }

    //  Store new deletable edges at first empty position or at the end of
    //  the array.
    if(insertPosition != bitsetsOfDeletableEdges->used) {
        insertArrayAtPos(bitsetsOfDeletableEdges, deletableEdges,
         insertPosition);
    }
    else {
        insertArray(bitsetsOfDeletableEdges, deletableEdges);
    }

    return 0;
}

//  Check if both of the other edges incident to x are not in deletableEdges. 
bool otherEdgesAreNonDeletable(bitset adjacencyList[], int numberOfVertices, 
 int x, int y, bitset deletableEdges, int edgeNumbering[][numberOfVertices]) {
    forEach(element, adjacencyList[x]) {
        if(element == y) {
            continue;
        }
        if(contains(deletableEdges, edgeNumbering[x][element])) {
            return false;
        }
    }
    return true;
}

// Add edges according to the three rules, return false if adding leads to
// contradiction
bool canAddNewArc(bitset adjacencyList[], int numberOfVertices,
 struct diGraph *orientation, int x, int y, bitset deletableEdges, 
 int edgeNumbering[][numberOfVertices]) {
    
    //  If the edge already exists, there cannot be any contradictions in the
    //  orientation
    if(contains(orientation->adjacencyList[x], y)) {
        return true;
    }

    if(contains(orientation->adjacencyList[y], x)) {
        return false;
    }

    if(size(orientation->adjacencyList[x]) >= 2) {
        return false;
    } 
    if(size(orientation->reverseAdjacencyList[y]) >= 2) {
        return false;
    } 

    //  Deletable edges incident to same vertex need to be one incoming, one
    //  outgoing.
    if(contains(deletableEdges, edgeNumbering[x][y])) {
        forEach(element, adjacencyList[x]) {
            if(element == y) {
                continue;
            }
            if(contains(deletableEdges, edgeNumbering[x][element])) {
                if(contains(orientation->adjacencyList[x], element)) {
                    return false;
                }
            }
        }
        forEach(element, adjacencyList[y]) {
            if(element == x) {
                continue;
            }
            if(contains(deletableEdges, edgeNumbering[y][element])) {
                if(contains(orientation->reverseAdjacencyList[y], element)) {
                    return false;
                }
            }
        }
    }
    else { // xy is not in deletableEdges

        // If xy was not deletable, it needs to be deletable in the current
        // orientation, i.e. x needs to have one incoming and one outgoing
        // apart from xy.
        if(size(orientation->adjacencyList[x]) >= 2 || 
         size(orientation->reverseAdjacencyList[x]) >= 2) {
            return false;
        }
        if(size(orientation->adjacencyList[y]) >= 2 ||
         size(orientation->reverseAdjacencyList[y]) >= 2) {
            return false;
        }

        //  If xy is not deletable then it needs the be oriented in opposite
        //  direction of other non-deletable edge incident to x.
        forEach(element, adjacencyList[x]) {
            if(element == y) {
                continue;
            }
            if(!contains(deletableEdges, edgeNumbering[x][element])) {
                if(contains(orientation->reverseAdjacencyList[x], y)) {
                    return false;
                }
                break;
            }
        }

        //  If xy is not deletable then it needs the be oriented in opposite
        //  direction of other non-deletable edge incident to y
        forEach(element, adjacencyList[y]) {
            if(element == x) {
                continue;
            }
            if(!contains(deletableEdges, edgeNumbering[y][element])) {
                if(contains(orientation->adjacencyList[y], x)) {
                    return false;
                }
                break;
            }
        }

    }
    addArc(orientation, x, y);

    //  If x has two outgoing and no incoming, add the final incoming.
    if(size(orientation->adjacencyList[x]) == 2 &&
     size(orientation->reverseAdjacencyList[x]) < 1) {
        int lastNeighbour = next(difference(adjacencyList[x],
         orientation->adjacencyList[x]), -1);
        if(!canAddNewArc(adjacencyList, numberOfVertices, orientation,
         lastNeighbour, x, deletableEdges, edgeNumbering)) {
            return false;
        }
    }

    //  If y has no outgoing and two incoming, add the final outgoing.
    if(size(orientation->adjacencyList[y]) == 0 &&
     size(orientation->reverseAdjacencyList[y]) == 2) {
        int lastNeighbour = next(difference(adjacencyList[y],
         orientation->reverseAdjacencyList[y]), -1);
        if(!canAddNewArc(adjacencyList, numberOfVertices, orientation, y,
         lastNeighbour, deletableEdges, edgeNumbering)) {
            return false;
        }
    }

    //  Deletable edges incident to same vertex need to be one incoming, one
    //  outgoing
    if(contains(deletableEdges, edgeNumbering[x][y])) {
        forEach(element, adjacencyList[x]) {
            if(element == y) {
                continue;
            }
            if(contains(deletableEdges, edgeNumbering[x][element])) {
                if(!canAddNewArc(adjacencyList, numberOfVertices, orientation,
                 element, x, deletableEdges, edgeNumbering)) {
                    return false;
                }
            }
        }
        forEach(element, adjacencyList[y]) {
            if(element == x) {
                continue;
            }
            if(contains(deletableEdges, edgeNumbering[y][element])) {
                if(!canAddNewArc(adjacencyList, numberOfVertices, orientation,
                 y, element, deletableEdges, edgeNumbering)) {
                    return false;
                }
            }
        }

        //  If one deletable edge and two nondeletable, the nondeletable need to
        //  be opposite of deletable.
        if(otherEdgesAreNonDeletable(adjacencyList, numberOfVertices, x, y,
         deletableEdges, edgeNumbering)) {
            forEach(element, adjacencyList[x]) {
                if(element == y) {
                    continue;
                }
                if(!canAddNewArc(adjacencyList, numberOfVertices, orientation,
                 element, x, deletableEdges, edgeNumbering)) {
                    return false;
                }
            }
        }

        if(otherEdgesAreNonDeletable(adjacencyList, numberOfVertices, y, x,
         deletableEdges, edgeNumbering)) {
            forEach(element, adjacencyList[y]) {
                if(element == x) {
                    continue;
                }
                if(!canAddNewArc(adjacencyList, numberOfVertices, orientation,
                 y, element, deletableEdges, edgeNumbering)) {
                    return false;
                }
            }
        }
    }
    else { // xy is not in deletableEdges

        // xy needs to be deletable, so if we have two incoming of y, we need an
        // outgoing.
        if(size(orientation->adjacencyList[y]) == 0 &&
         size(orientation->reverseAdjacencyList[y]) == 2) {
            int lastNeighbour = next(difference(adjacencyList[y],
             orientation->adjacencyList[y]), -1);
            if(!canAddNewArc(adjacencyList, numberOfVertices, orientation, y,
             lastNeighbour, deletableEdges, edgeNumbering)) {
                return false;
            }
        }

        // xy needs to be deletable, so if we have one outgoing, one incoming to
        // y, we need an incoming.
        if(size(orientation->adjacencyList[y]) == 1 &&
         size(orientation->reverseAdjacencyList[y]) == 1) {
            int lastNeighbour = next(difference(adjacencyList[y],
             union(orientation->adjacencyList[y], 
             orientation->reverseAdjacencyList[y])), -1);
            if(!canAddNewArc(adjacencyList, numberOfVertices, orientation,
             lastNeighbour, y, deletableEdges, edgeNumbering)) {
                return false;
            }
        }

        //  If xy is not deletable then it needs the be oriented in opposite
        //  direction of other non-deletable edge incident to x
        forEach(element, adjacencyList[x]) {
            if(element == y) {continue; } if(!contains(deletableEdges, 
             edgeNumbering[x][element])) {
                if(!canAddNewArc(adjacencyList, numberOfVertices, orientation,
                 x, element, deletableEdges, edgeNumbering)) {
                    return false;
                }
                break;
            }
        }

        //  If xy is not deletable then it needs the be oriented in opposite
        //  direction of other non-deletable edge incident to y
        forEach(element, adjacencyList[y]) {
            if(element == x) {
                continue;
            }
            if(!contains(deletableEdges, edgeNumbering[y][element])) {
                if(!canAddNewArc(adjacencyList, numberOfVertices, orientation,
                 element, y, deletableEdges, edgeNumbering)) {
                    return false;
                }
                break;
            }
        }

    }
    return true;
}

//  Loop over all edges and try orienting them in both directions.
bool canCompleteCompOrientation(bitset adjacencyList[], int numberOfVertices,
 struct options *options, struct diGraph *orientation, bitset deletableEdges,
 int edgeNumbering[][numberOfVertices], int endpoint1, int endpoint2) {

    if(endpoint2 == -1 && endpoint1 < (numberOfVertices - 1)) {
        return canCompleteCompOrientation(adjacencyList, numberOfVertices,
         options, orientation, deletableEdges, edgeNumbering, endpoint1 + 1,
         next(adjacencyList[endpoint1 + 1], endpoint1 + 1));
    }

    //  We have oriented all edges.
    if(endpoint2 == -1 && endpoint1 == numberOfVertices - 1) {
        if(orientation->numberOfArcs != 3*numberOfVertices/2) {
            fprintf(stderr, "%s\n", "Something went wrong");
        }

        //  Check if formed orientation actually is complementary.
        bitset complementDeletableEdges = getDeletableEdges(orientation,
         numberOfVertices, edgeNumbering);
        if(equals(union(deletableEdges, complementDeletableEdges),
         complement(EMPTY, 3*numberOfVertices/2))) {
            if(options->printFlag) {
                printDeletableEdges(numberOfVertices, edgeNumbering,
                 orientation->adjacencyList, complementDeletableEdges);
                printDiGraph(orientation);
            }
            return true;
        }
        return false;
    }

    //  If already oriented, go to next edge.
    if(contains(orientation->adjacencyList[endpoint1], endpoint2) ||
     contains(orientation->adjacencyList[endpoint2], endpoint1)) {
        return canCompleteCompOrientation(adjacencyList, numberOfVertices,
         options, orientation, deletableEdges, edgeNumbering, endpoint1,
         next(adjacencyList[endpoint1], endpoint2));
    }

    //  Make copy of orientation
    struct diGraph orientationCopy = {.numberOfVertices = numberOfVertices};
    orientationCopy.adjacencyList = malloc(sizeof(bitset)*numberOfVertices);
    orientationCopy.reverseAdjacencyList = 
     malloc(sizeof(bitset)*numberOfVertices);
    memcpy(orientationCopy.adjacencyList, orientation->adjacencyList,
     sizeof(bitset)*numberOfVertices);
    memcpy(orientationCopy.reverseAdjacencyList,
     orientation->reverseAdjacencyList, sizeof(bitset)*numberOfVertices);
    orientationCopy.numberOfArcs = orientation->numberOfArcs;

    //  Try adding endpoint1->endpoint2
    if(canAddNewArc(adjacencyList, numberOfVertices, orientation, endpoint1,
     endpoint2, deletableEdges, edgeNumbering)) {

        //  Continue with next edge.
        if(canCompleteCompOrientation(adjacencyList, numberOfVertices, options,
         orientation, deletableEdges, edgeNumbering, endpoint1, 
         next(adjacencyList[endpoint1], endpoint2))) {
            free(orientationCopy.adjacencyList);
            free(orientationCopy.reverseAdjacencyList);
            return true;
       }
    }

    //  Put orientation back to before we added endpoint1->endpoint2.
    memcpy(orientation->adjacencyList, orientationCopy.adjacencyList,
     sizeof(bitset)*numberOfVertices);
    memcpy(orientation->reverseAdjacencyList,
     orientationCopy.reverseAdjacencyList, sizeof(bitset)*numberOfVertices);
    orientation->numberOfArcs = orientationCopy.numberOfArcs;

    //  Try adding endpoint2->endpoint1.
    if(canAddNewArc(adjacencyList, numberOfVertices, orientation, endpoint2,
     endpoint1, deletableEdges, edgeNumbering)) {

        //  Continue with next edge.
        if(canCompleteCompOrientation(adjacencyList, numberOfVertices, options,
         orientation, deletableEdges, edgeNumbering, endpoint1,
         next(adjacencyList[endpoint1], endpoint2))) {
            free(orientationCopy.adjacencyList);
            free(orientationCopy.reverseAdjacencyList);
            return true;
        }
    }

    //  Both orientations lead to contradiction.
    free(orientationCopy.adjacencyList);
    free(orientationCopy.reverseAdjacencyList);
    return false;
}

bool hasComplementaryOrientation(bitset adjacencyList[], int numberOfVertices,
 struct options *options, bitset deletableEdgesOfOrientationTocomplement, 
 int edgeNumbering[][numberOfVertices]) {

    //  This will complement the given orientation.
    struct diGraph orientation = {.numberOfVertices = numberOfVertices};
    orientation.adjacencyList = malloc(sizeof(bitset)*numberOfVertices);
    orientation.reverseAdjacencyList = malloc(sizeof(bitset)*numberOfVertices);
    emptyGraph(&orientation);

    //  Fix a first arc, does not matter which or in what direction.
    //  (Orientation with opposite arcs has same deletable edges.)
    if(!canAddNewArc(adjacencyList, numberOfVertices, &orientation, 0,
     next(adjacencyList[0], -1), deletableEdgesOfOrientationTocomplement,
      edgeNumbering)) {
        return false;
    }

    bool hasCompOrientation = canCompleteCompOrientation(adjacencyList,
     numberOfVertices, options, &orientation,
     deletableEdgesOfOrientationTocomplement, edgeNumbering, 0,
     next(adjacencyList[0], -1));

    free(orientation.adjacencyList);
    free(orientation.reverseAdjacencyList);
    return hasCompOrientation;
}

//  Generate strong orientations of graph, get deletable edges and perform one
//  of the exact methods.
int generateAllOrientations(bitset adjacencyList[], struct options *options,
 struct counters *numberOf, int numberOfVertices, 
 int edgeNumbering[][numberOfVertices], Array *bitsetsOfDeletableEdges, 
 struct diGraph *orientation, int endpoint1, int endpoint2) {

    int frankNumberUpperBound = 0;
    if(endpoint2 == -1 && endpoint1 < (numberOfVertices - 1)) {
        frankNumberUpperBound = generateAllOrientations(adjacencyList, options,
         numberOf, numberOfVertices, edgeNumbering, bitsetsOfDeletableEdges,
         orientation, endpoint1 + 1, next(adjacencyList[endpoint1 + 1], 
         endpoint1 + 1));
        return frankNumberUpperBound;
    }

    //  All edges are oriented.
    if(endpoint2 == -1 && endpoint1 == numberOfVertices - 1) {

        numberOf->totalOrientationsGenerated++;

        //  Skip orientations which do not have correct remainder.
        if(options->singleGraphFlag) {
            if(numberOf->totalOrientationsGenerated % options->modulo != 
             options->remainder) {
               return 0; 
            }
        }

        if(!isStronglyConnected(orientation)) {
            return 0;
        }

        bitset deletableEdges = getDeletableEdges(orientation, numberOfVertices,
         edgeNumbering);

        //  Check if there is a vertex with three non-deletable incident edges.
        //  In this case orientation has no complementary orientation giving
        //  fn=2.
        for(int i = 0; i < numberOfVertices; i++) {
            bool noIncidentEdgesDeletable = true;
            forEach(nbr, adjacencyList[i])  {
                if(contains(deletableEdges, edgeNumbering[i][nbr])) {
                    noIncidentEdgesDeletable = false;
                }
            }
            if(noIncidentEdgesDeletable) {
                return 0;
            }
        }

        numberOf->generatedOrientations++;

        //  Try finding a complement to the current orientation.
        if(!options->bruteForceFlag) {
            if(hasComplementaryOrientation(adjacencyList, numberOfVertices,
             options, deletableEdges, edgeNumbering)) {
                if(options->printFlag) {
                    printDeletableEdges(numberOfVertices, edgeNumbering,
                     orientation->adjacencyList, deletableEdges);
                    printDiGraph(orientation);
                }
                return 2;
            } 
            return 0;
        }

        //  If not complementFlag, try using the bruteforce method of comparing
        //  all orientations pairwise.
        return getIntermediateFrankNumber(options, numberOf, numberOfVertices,
         edgeNumbering, bitsetsOfDeletableEdges, deletableEdges);
    }

    //  Orient edge and continue with next edge.
    addArc(orientation, endpoint1, endpoint2);
    if(size(orientation->adjacencyList[endpoint1]) != 3 &&
     size(orientation->reverseAdjacencyList[endpoint2]) != 3) {
        frankNumberUpperBound = generateAllOrientations(adjacencyList, options,
         numberOf, numberOfVertices, edgeNumbering, bitsetsOfDeletableEdges, 
         orientation, endpoint1, next(adjacencyList[endpoint1], endpoint2));
    }
    removeArc(orientation, endpoint1, endpoint2);

    if(frankNumberUpperBound) {
        return frankNumberUpperBound;
    }

    //  Orient edge in other way and continue.
    addArc(orientation, endpoint2, endpoint1);
    if(size(orientation->reverseAdjacencyList[endpoint1]) != 3 && 
     size(orientation->adjacencyList[endpoint2]) != 3) {
        frankNumberUpperBound = generateAllOrientations(adjacencyList, options,
         numberOf, numberOfVertices, edgeNumbering, bitsetsOfDeletableEdges,
         orientation, endpoint1, next(adjacencyList[endpoint1], endpoint2));
    }
    removeArc(orientation, endpoint2, endpoint1);

    if(frankNumberUpperBound) {
        return frankNumberUpperBound;
    }

    //  None of the orientations of this edge led to an orientation of the graph
    //  which has a second orientation giving fn=2.
    return 0;
}

int findFrankNumber(bitset adjacencyList[], int numberOfVertices,
 struct options *options, struct counters *numberOf) {
    Array bitsetsOfDeletableEdges;
    initArray(&bitsetsOfDeletableEdges, options->sizeOfArray);

    int edgeNumbering[numberOfVertices][numberOfVertices];
    numberEdges(adjacencyList, numberOfVertices, edgeNumbering);

    struct diGraph orientation = {.numberOfVertices = numberOfVertices};
    orientation.adjacencyList = malloc(sizeof(bitset)*numberOfVertices);
    orientation.reverseAdjacencyList = malloc(sizeof(bitset)*numberOfVertices);
    emptyGraph(&orientation);

    int frankNumber = generateAllOrientations(adjacencyList, options, numberOf,
     numberOfVertices, edgeNumbering, &bitsetsOfDeletableEdges, &orientation,
     -1, -1);

    //  In bruteforce case, we now have a list of bitsets corresponding to
    //  deletable edges of (all) orientations.
    if(options->bruteForceFlag) {
        numberOf->storedBitsets = bitsetsOfDeletableEdges.used;
        if(numberOf->storedBitsets > options->sizeOfArray) {
            options->sizeOfArray = bitsetsOfDeletableEdges.size;
        }
        if(options->verboseFlag) {
            fprintf(stderr, "\tBitsets stored: %llu, size of array %llu\n", 
             numberOf->storedBitsets, options->sizeOfArray);
        }

        //  Count empty bitsets stored and check that there are enough
        //  orientations for the Frank number to make sense. (This should of
        //  course always be the case.)
        bitset universe = EMPTY;
        for(size_t i = 0; i < bitsetsOfDeletableEdges.used; i++) {
            if(isEmpty(bitsetsOfDeletableEdges.array[i])) {
                numberOf->emptyBitsetsStored++;
            }
            universe = union(universe, bitsetsOfDeletableEdges.array[i]);
        }
        if(options->verboseFlag) {
            fprintf(stderr, "\tEmpty bitsets stored: %llu \n", 
             numberOf->emptyBitsetsStored);
        }
        if(!equals(universe, complement(EMPTY, 3*numberOfVertices/2))) {
            fprintf(stderr, "%s\n",
             "Error: Not enough orientations for Frank number to make sense.");
        }
    }

    freeArray(&bitsetsOfDeletableEdges);
    free(orientation.adjacencyList);
    free(orientation.reverseAdjacencyList);
    return frankNumber;
}

//******************************************************************************
//
//                              Heuristic algorithm
//
//******************************************************************************

//  Array and bitset representation of a cycle.
struct cycle {
    bitset cycleElements;
    int numberOfElements;
    int *cycle;
};

//  Count the odd cycles in a complement of the perfect F F. Assuming
//  graphs to be cubic and bridgeless. We also store for each even cycle a
//  maximal F in M.
bool containsTwoOddCycles(bitset adjacencyList[], int numberOfVertices,
 int F[], struct cycle oddCycles[], int M[]) {

    for(int i = 0; i < numberOfVertices; i++) {
        M[i] = -1;
    }
    int numberOfOddCycles = 0;
    bitset uncheckedVertices = complement(EMPTY, numberOfVertices);

    //  Loop over all cycles and check parity
    //  Store the odd edges of each cycle in M
    forEach(element, uncheckedVertices) {
        int currentVertex = element;
        int previousVertex = -1;
        bool cycleIsOdd = false;
        bitset cycle = EMPTY;
        if(numberOfOddCycles < 2) {
            oddCycles[numberOfOddCycles].numberOfElements = 0;
        }
        do {
            removeElement(uncheckedVertices, currentVertex);
            add(cycle, currentVertex);
            if(numberOfOddCycles < 2) {
                oddCycles[numberOfOddCycles].cycle[
                 oddCycles[numberOfOddCycles].numberOfElements++] = 
                 currentVertex;
            }
            int nextVertex = next(adjacencyList[currentVertex], -1);
            while(nextVertex == previousVertex ||
             nextVertex == F[currentVertex]) {
                nextVertex = next(adjacencyList[currentVertex], nextVertex);
            }
            if(M[currentVertex] == -1) {
                M[currentVertex] = nextVertex;
                M[nextVertex] = currentVertex;
            }
            previousVertex = currentVertex;
            currentVertex = nextVertex;
            cycleIsOdd = !cycleIsOdd;
        } while(currentVertex != element);

        if(cycleIsOdd) {
            if(numberOfOddCycles < 2) {
                oddCycles[numberOfOddCycles].cycleElements = cycle;
            }
            numberOfOddCycles++;
            if(numberOfOddCycles > 2) {
                return false;
            }
        }   
    }
    return numberOfOddCycles == 2;
}

//  Add maximal F of odd cycles - x1 - x2 to M.
void getOddCycleMatching(bitset adjacencyList[], int numberOfVertices,
 struct cycle oddCycles[], int indexOfx1, int indexOfx2, int M[]) {

    int currentIndex = indexOfx1;
    bool addToMatching = false;
    do { 
        int nextIndex = (currentIndex + 1) % oddCycles[0].numberOfElements;
        if(addToMatching) {
            M[oddCycles[0].cycle[nextIndex]] = oddCycles[0].cycle[currentIndex];
            M[oddCycles[0].cycle[currentIndex]] = oddCycles[0].cycle[nextIndex];
        }
        addToMatching = !addToMatching;
        currentIndex = nextIndex;
    } while(currentIndex != indexOfx1);

    currentIndex = indexOfx2;
    addToMatching = false;
    do { 
        int nextIndex = (currentIndex + 1) % oddCycles[1].numberOfElements;
        if(addToMatching) {
            M[oddCycles[1].cycle[nextIndex]] = oddCycles[1].cycle[currentIndex];
            M[oddCycles[1].cycle[currentIndex]] = oddCycles[1].cycle[nextIndex];
        }
        addToMatching = !addToMatching;
        currentIndex = nextIndex;
    } while(currentIndex != indexOfx2);
}

//  Find the index of u in arr[].
int findInArray(int u, int arr[], int arrLength) {
    for(int i = 0; i < arrLength; i++) {
        if(arr[i] == u) {
            return i;
        }
    }
    return -1;
}

// Check if orientation of F - {x1,x2,(y1,y2)} is consistent on the cycle
// containing u and v.
bool circuitOrientationIsConsistent(bitset adjacencyList[],
 int numberOfVertices, int M[], int F[], 
 int circuitOrientation[], int u, int v) {

    //  If circuit containing u of F - {x1,x2,(y1,y2)} not yet oriented, orient
    //  it.
    if(circuitOrientation[u] == -1) {

        //  Orient the edges incident to u such that they are consistent with
        //  the edges incident to v on the cycle containing u and v. If
        //  circuitOrientation[v] is still -1, the direction we orient the
        //  edges incident with u does not matter. 
        int takeMaximalMatching = (circuitOrientation[v] == F[v]);
        int currentVertex = u;
        do {
            int nextVertex = takeMaximalMatching ? M[currentVertex] : 
             F[currentVertex];
            circuitOrientation[currentVertex] = nextVertex;
            currentVertex = nextVertex;
            takeMaximalMatching = !takeMaximalMatching;
        } while (currentVertex != u);
    }

    if(circuitOrientation[v] == -1) {
        int takeMaximalMatching = (circuitOrientation[u] == F[u]);
        int currentVertex = v;
        do {
            int nextVertex = takeMaximalMatching ? M[currentVertex] :
             F[currentVertex];
            circuitOrientation[currentVertex] = nextVertex;
            currentVertex = nextVertex;
            takeMaximalMatching = !takeMaximalMatching;
        } while (currentVertex != v);
    }
    return (circuitOrientation[u] == F[u]) == (circuitOrientation[v] == M[v]);
}

//  If we are checking the heuristic with the even cycle, M might not be a
//  correct maximal matching of this even cycle. Redo it.
void rematch(bitset adjacencyList[], int numberOfVertices, int M[], int F[],
 int y1, int y2) {
    int previousVertex = y2;
    int currentVertex = y1;
    bool addToMaximalMatching = false;
    do {
        int nextVertex = next(difference(adjacencyList[currentVertex], 
         union(singleton(F[currentVertex]), singleton(previousVertex))), -1);
        if(addToMaximalMatching) {
            M[currentVertex] =  nextVertex;
            M[nextVertex] = currentVertex;
        }
        previousVertex = currentVertex;
        currentVertex = nextVertex;
        addToMaximalMatching = !addToMaximalMatching; 
    } while(currentVertex != y2);

    M[y1] = y2;
    M[y2] = y1;
}

//  For checking cyclic connectivity.
void DFS(bitset adjacencyList[], int numberOfVertices, bitset *component,
 bitset *uncheckedVertices, int v, int parent, bool *cycleFound) {

    //  If checked before: cycle found.
    if(contains(*component, v)) {
        *cycleFound = true;
        return;
    }
    removeElement(*uncheckedVertices, v);
    add(*component, v);

    //  Do not go back to parent.
    forEach(nbr, difference(adjacencyList[v], singleton(parent))) {
        DFS(adjacencyList, numberOfVertices, component, uncheckedVertices, nbr,
         v, cycleFound);
    } 
}

bool isCyclicallyConnected(bitset adjacencyList[], int numberOfVertices) {
    bitset uncheckedVertices = complement(EMPTY, numberOfVertices);
    int numberOfComponents = 0;
    int numberOfComponentsWithCycle = 0;
    bitset components[numberOfVertices];
    forEach(v, uncheckedVertices) {
        numberOfComponents++;
        components[numberOfComponents - 1] = EMPTY;
        bool cycleFound = false;
        DFS(adjacencyList, numberOfVertices, &components[numberOfComponents - 1],
         &uncheckedVertices, v, -1, &cycleFound);
        if(cycleFound) {
            numberOfComponentsWithCycle++;
        }
        if(numberOfComponentsWithCycle >= 2) {
            return false;
        }
    }
    return true;
}

#define removeEdgeFromAdjList(adjacencyList, endpoint1, endpoint2) {\
    removeElement(adjacencyList[endpoint1], endpoint2);\
    removeElement(adjacencyList[endpoint2],endpoint1);\
}
#define addEdgeToAdjList(adjacencyList, endpoint1, endpoint2) {\
    add(adjacencyList[endpoint1], endpoint2);\
    add(adjacencyList[endpoint2],endpoint1);\
}

//  Check if edge is a strong 2-edge under the assumption it is valuated 2 in
//  the flow. Hence, we only check it is not part of some cycle-separating
//  3-edge-set containing two other edges from circuitOrientation (This is a
//  sufficient condition.)
bool edgeIsStrong2Edge(bitset adjacencyList[], int numberOfVertices,
 int endpoint1, int endpoint2, int circuitOrientation[]) {
    bool hasCyclic211cut = false;
    removeEdgeFromAdjList(adjacencyList, endpoint1, endpoint2);

    //  Loop over all pairs of edges in the perfect F.
    for(int i = 0; i < numberOfVertices; i++) {

        //  Should not look at edge we already suppressed.
        if(circuitOrientation[i] == -1) {
            continue;
        }
        removeEdgeFromAdjList(adjacencyList, i, circuitOrientation[i]);

        for(int j = i+1; j < numberOfVertices; j++) {
            if(circuitOrientation[j] == -1) {
                continue;
            }
            removeEdgeFromAdjList(adjacencyList, j, circuitOrientation[j]);

            if(!isCyclicallyConnected(adjacencyList, numberOfVertices)) {
                hasCyclic211cut = true;
            }
            addEdgeToAdjList(adjacencyList, j, circuitOrientation[j]);

            if(hasCyclic211cut) {
                break;
            }
        }

        addEdgeToAdjList(adjacencyList, i, circuitOrientation[i]);
        if(hasCyclic211cut) {
            break;
        }
    }
    addEdgeToAdjList(adjacencyList, endpoint1, endpoint2);
    return !hasCyclic211cut;
}

//  Are the suppressed strong 2-edges in the nz 4-flow deletable?
bool suppressedEdgesAreDeletable(bitset adjacencyList[], int numberOfVertices,
 int circuitOrientation[], int edgesBetweenCycles[],
 int numberOfEdgesBetweenCycles) {
    bool edgesAreDeletable = true;
    for(int i = 0; i < numberOfEdgesBetweenCycles; i++) {
        removeEdgeFromAdjList(adjacencyList, edgesBetweenCycles[2*i],
         edgesBetweenCycles[2*i+1]);
    }
    for(int i = 0; i < numberOfEdgesBetweenCycles; i++) {
        if(!edgeIsStrong2Edge(adjacencyList, numberOfVertices,
         edgesBetweenCycles[2*i], next(adjacencyList[edgesBetweenCycles[2*i]], 
         -1), circuitOrientation)){
            edgesAreDeletable = false;
            break;
        }
        if(!edgeIsStrong2Edge(adjacencyList, numberOfVertices,
         edgesBetweenCycles[2*i+1], next(
         adjacencyList[edgesBetweenCycles[2*i+1]], -1), circuitOrientation)){
            edgesAreDeletable = false;
            break;
        }
    }
    for(int i = 0; i < numberOfEdgesBetweenCycles; i++) {
        addEdgeToAdjList(adjacencyList, edgesBetweenCycles[2*i],
         edgesBetweenCycles[2*i+1]);
    }
    return edgesAreDeletable;
}

//  Used for double checking heuristic algorithm.
void verifyOddnessHeuristicOrientations(bitset adjacencyList[],
 int numberOfVertices, struct options *options, int circuitOrientation[], 
 int F[], int M[], int edgesBetweenCycles[], int numberOfEdgesBetweenCycles); 

// Generate all perfect matchings of the graph and check for each of the
// complementary 2-factors whether one of the configurations for the sufficient
// conditions are present.
bool hasSufficientCondition(bitset adjacencyList[], int numberOfVertices,
 struct options *options, struct counters *numberOf, bitset remainingVertices,
 int F[]) {

    //  If this holds, F is a perfect matching.
    int nextVertex = next(remainingVertices, -1);
    if(nextVertex == -1) {

        struct cycle oddCycles[2];
        oddCycles[0].cycle = malloc(sizeof(int)*numberOfVertices);
        oddCycles[1].cycle = malloc(sizeof(int)*numberOfVertices);
        int M[numberOfVertices];
        if(containsTwoOddCycles(adjacencyList, numberOfVertices, F, oddCycles,
         M)) {

            //  Check if odd cycles are connected via an edge uv
            forEach(u, oddCycles[0].cycleElements) {
                int v = F[u];
                if(contains(oddCycles[1].cycleElements, v)) {
                    int indexOfx1 = findInArray(u, oddCycles[0].cycle,
                     oddCycles[0].numberOfElements);
                    int indexOfx2 = findInArray(v, oddCycles[1].cycle,
                     oddCycles[1].numberOfElements);

                    //  Add a maximal matching of the odd cycles to the maximal
                    //  matching M of G - F
                    getOddCycleMatching(adjacencyList, numberOfVertices,
                     oddCycles, indexOfx1, indexOfx2, M);

                    int u1 = oddCycles[0].cycle[
                        (indexOfx1 + 1) % oddCycles[0].numberOfElements];
                    int u2 = oddCycles[1].cycle[
                        (indexOfx2 + 1) % oddCycles[1].numberOfElements];
                    int v1 = oddCycles[0].cycle[
                        (oddCycles[0].numberOfElements + indexOfx1 - 1) % 
                        oddCycles[0].numberOfElements];
                    int v2 = oddCycles[1].cycle[
                        (oddCycles[1].numberOfElements + indexOfx2 - 1) % 
                        oddCycles[1].numberOfElements];                    
                    
                    //  Orient cycles of F and check condition
                    int circuitOrientation[numberOfVertices];
                    //  Can be optimized!
                    for(int i = 0; i < numberOfVertices; i++) {
                        circuitOrientation[i] = -1;
                    }
                    if(circuitOrientationIsConsistent(adjacencyList,
                     numberOfVertices, M, F, circuitOrientation,
                     u1, v1) && 
                     circuitOrientationIsConsistent(adjacencyList, 
                     numberOfVertices, M, F, circuitOrientation, 
                     u2, v2)) {
                        int edgesBetweenCycles[] = {u,v};
                        if(suppressedEdgesAreDeletable(adjacencyList,
                         numberOfVertices, circuitOrientation,
                         edgesBetweenCycles, 1)) {
                            numberOf->graphsSatisfyingFirstOddness++;
                            if(options->doublecheckFlag || options->printFlag) {
                                verifyOddnessHeuristicOrientations(
                                 adjacencyList, numberOfVertices, options,
                                 circuitOrientation, F, M, edgesBetweenCycles,
                                 1);
                            }
                            free(oddCycles[0].cycle);
                            free(oddCycles[1].cycle);
                            return true;
                        }
                        if(options->verboseFlag) {
                            fprintf(stderr, "Not deletable: first\n");
                        }
                    }
                    continue;
                }
                if(!contains(oddCycles[0].cycleElements, v)) {
                    int nbrOfU = v;
                    forEach(nbrOfV, adjacencyList[nbrOfU]) {
                        if(nbrOfV == u) {
                            continue;
                        }
                        v = next(intersection(adjacencyList[nbrOfV],
                         oddCycles[1].cycleElements),-1);
                        if(v == -1) {
                            continue;
                        }
                        int indexOfx1 = findInArray(u, oddCycles[0].cycle,
                         oddCycles[0].numberOfElements);
                        int indexOfx2 = findInArray(v, oddCycles[1].cycle,
                         oddCycles[1].numberOfElements);
                        getOddCycleMatching(adjacencyList, numberOfVertices,
                         oddCycles, indexOfx1, indexOfx2, M);
                        int u1 = oddCycles[0].cycle[(indexOfx1 + 1) %
                         oddCycles[0].numberOfElements];
                        int u2 = oddCycles[1].cycle[(indexOfx2 + 1) %
                         oddCycles[1].numberOfElements];
                        int v1 = oddCycles[0].cycle[
                         (oddCycles[0].numberOfElements + indexOfx1 - 1) % 
                         oddCycles[0].numberOfElements];
                        int v2 = oddCycles[1].cycle[
                         (oddCycles[1].numberOfElements + indexOfx2 - 1) %
                         oddCycles[1].numberOfElements];
                        int w1 = next(difference(adjacencyList[nbrOfU],
                         union(singleton(nbrOfV), singleton(F[nbrOfU]))),-1);  
                        int w2 = next(difference(adjacencyList[nbrOfV],
                         union(singleton(nbrOfU), singleton(F[nbrOfV]))), -1);     
                        
                        //  Orient cycles and check condition
                        int circuitOrientation[numberOfVertices];
                        for(int i = 0; i < numberOfVertices; i++) {
                            circuitOrientation[i] = -1;
                        }

                        //  Adapt the matching of the even cycle such that M is
                        //  still maximal in C - {x1,x2,y1,y2}
                        if(M[nbrOfU] != nbrOfV) {
                            rematch(adjacencyList, numberOfVertices, M, F,
                             nbrOfU, nbrOfV);
                        }

                        //  Check if orientations are consistent
                        if(circuitOrientationIsConsistent(adjacencyList,
                         numberOfVertices, M, F, circuitOrientation, u1, v1) && 
                         circuitOrientationIsConsistent(adjacencyList, 
                         numberOfVertices, M, F, circuitOrientation, u2, v2) && 
                         circuitOrientationIsConsistent(adjacencyList, 
                         numberOfVertices, M, F, circuitOrientation, w1, w2)) {
                            int edgesBetweenCycles[] = {u, nbrOfU, nbrOfV, v};
                            if(suppressedEdgesAreDeletable(adjacencyList,
                             numberOfVertices, circuitOrientation,
                             edgesBetweenCycles, 2)) {
                                numberOf->graphsSatisfyingSecondOddness++;
                                if(options->doublecheckFlag || 
                                 options->printFlag) {
                                    verifyOddnessHeuristicOrientations(
                                     adjacencyList, numberOfVertices, options,
                                     circuitOrientation, F, M,
                                     edgesBetweenCycles, 2);
                                }
                                free(oddCycles[0].cycle);
                                free(oddCycles[1].cycle);
                                return true;
                            }
                            if(options->verboseFlag) {
                                fprintf(stderr, "Not deletable\n");
                            }
                        }
                        continue;
                    }
                }
            }
        }

        //  None of the configurations were present for this perfect matching.
        free(oddCycles[0].cycle);
        free(oddCycles[1].cycle);
        return false;
    }

    //  F is not yet a perfect matching here. 
    forEach(neighbor, intersection(adjacencyList[nextVertex],
     remainingVertices)) {
        F[neighbor] = nextVertex;
        F[nextVertex] = neighbor;
        bitset newRemainingVertices = difference(remainingVertices,
         union(singleton(nextVertex), singleton(neighbor)));
        if(hasSufficientCondition(adjacencyList, numberOfVertices, options,
         numberOf, newRemainingVertices, F)) {
            return true;
        }
    }
    return false;
}

//  Make the concrete orientations for double checking the heuristic algorithm.
void orient2FactorCyclesInComplementaryOrientations(bitset adjacencyList[],
 int F[], int circuitOrientation[], int startingVertex,
 bitset *uncheckedVertices, struct diGraph *orientation1,
 struct diGraph *orientation2) {
    int currentVertex = startingVertex;

    //  Currentvertex lies on a cycle of the 2-factor. Circuitorientation
    //  orients an outgoing edge of this cycle from one of the neighbours of
    //  currentvertex. We want to orient the cycle in this direction. Hence,
    //  previousvertex should be oriented in the circuitorientation and
    //  prev->curr should be the direction we are orienting, hence
    //  circuitOrientation[prev] should be F[prev];
    int previousVertex = next(difference(adjacencyList[currentVertex],
     singleton(F[currentVertex])),-1);
    if(circuitOrientation[previousVertex] == -1 || 
     circuitOrientation[previousVertex] != F[previousVertex]) {
        previousVertex = next(difference(adjacencyList[currentVertex],
         singleton(F[currentVertex])), previousVertex);
    }
    do {
        removeElement((*uncheckedVertices), currentVertex);
        int nextVertex = next(adjacencyList[currentVertex], -1);
        while(nextVertex == previousVertex || nextVertex == F[currentVertex]) {
            nextVertex = next(adjacencyList[currentVertex], nextVertex);
        }
        if(circuitOrientation[nextVertex] == currentVertex) {
            addArc(orientation2, currentVertex, nextVertex);
            removeArc(orientation2, nextVertex, currentVertex);
        }
        else if(circuitOrientation[currentVertex] != nextVertex && 
         circuitOrientation[nextVertex] != currentVertex) {
            addArc(orientation1, currentVertex, nextVertex);
            addArc(orientation2, currentVertex, nextVertex);
        }
        previousVertex = currentVertex;
        currentVertex = nextVertex;
    } while(currentVertex != startingVertex);
}

// Make the concrete orientations for double checking the heuristic algorithm.
void verifyOddnessHeuristicOrientations(bitset adjacencyList[],
 int numberOfVertices, struct options *options, int circuitOrientation[],
 int F[], int M[], int edgesBetweenCycles[], int numberOfEdgesBetweenCycles) {

    struct diGraph orientation1 = {.numberOfVertices = numberOfVertices, 
     .numberOfArcs = 0};
    orientation1.adjacencyList = malloc(sizeof(bitset)*numberOfVertices);
    orientation1.reverseAdjacencyList = malloc(sizeof(bitset)*numberOfVertices);
    emptyGraph(&orientation1);
    struct diGraph orientation2 = {.numberOfVertices = numberOfVertices, 
     .numberOfArcs = 0};
    orientation2.adjacencyList = malloc(sizeof(bitset)*numberOfVertices);
    orientation2.reverseAdjacencyList = malloc(sizeof(bitset)*numberOfVertices);
    emptyGraph(&orientation2);

    // Add arc between u and v and add endpoints to bitset.
    bitset endpoints = EMPTY;
    for(int i = 0; i < numberOfEdgesBetweenCycles; i++) {
        addArc(&orientation1, edgesBetweenCycles[2*i], 
         edgesBetweenCycles[2*i+1]);
        addArc(&orientation2, edgesBetweenCycles[2*i+1], 
         edgesBetweenCycles[2*i]);
        add(endpoints, edgesBetweenCycles[2*i]);
        add(endpoints, edgesBetweenCycles[2*i+1]);
    }

    //  Add arcs from the circuitOrientation
    for(int i = 0; i < numberOfVertices; i++) {

        //  If i is one of the endpoints it does not belong to the circuits.
        if(contains(endpoints, i)) {
            continue;
        }

        //  Some part of the circuits might not yet be oriented. Do this now.
        if(circuitOrientation[i] == -1) {
            int takeMaximalMatching = true;
            int currentVertex = i;
            do {
                // fprintf(stderr, "%d\n", currentVertex);
                int nextVertex = takeMaximalMatching ? M[currentVertex] :
                 F[currentVertex];
                circuitOrientation[currentVertex] = nextVertex;
                currentVertex = nextVertex;
                takeMaximalMatching = !takeMaximalMatching;
            } while (currentVertex != i);
        }
        addArc(&orientation1, circuitOrientation[i], i);
        addArc(&orientation2, i, circuitOrientation[i]);
    }

    //  Orient 2-factor cycles
    bitset uncheckedVertices = complement(EMPTY, numberOfVertices);

    // Start orienting cycle at every endpoint of edge on cycle.
    for(int i = 0; i < 2*numberOfEdgesBetweenCycles; i++) { 
        if(contains(uncheckedVertices, edgesBetweenCycles[i])) {
            orient2FactorCyclesInComplementaryOrientations(adjacencyList, F, 
             circuitOrientation, edgesBetweenCycles[i], &uncheckedVertices,
             &orientation1, &orientation2);
        }
    }
    forEach(element, uncheckedVertices) {
        orient2FactorCyclesInComplementaryOrientations(adjacencyList, F,
         circuitOrientation, element, &uncheckedVertices, &orientation1,
         &orientation2);
    }

    if(!isStronglyConnected(&orientation1) || 
     !isStronglyConnected(&orientation2)) {
        fprintf(stderr, 
         "Error: orientations from oddness 2 heuristic not strongly connected!\n");
        exit(1);
    }

    int edgeNumbering[numberOfVertices][numberOfVertices];
    numberEdges(adjacencyList, numberOfVertices, edgeNumbering);
    bitset deletableEdges1 = getDeletableEdges(&orientation1, numberOfVertices, 
     edgeNumbering);
    bitset deletableEdges2 = getDeletableEdges(&orientation2, numberOfVertices, 
     edgeNumbering);

    if(options->printFlag) {
        printDeletableEdges(numberOfVertices, edgeNumbering,
         orientation1.adjacencyList, deletableEdges1);
        printDiGraph(&orientation1);
        printDeletableEdges(numberOfVertices, edgeNumbering, 
         orientation2.adjacencyList, deletableEdges2);
        printDiGraph(&orientation2);
    }

    if(!(equals(union(deletableEdges1, deletableEdges2),
     complement(EMPTY, 3*numberOfVertices/2)))) {
        fprintf(stderr, 
         "Error: orientations from oddness 2 heuristic are not complementary!\n");
        exit(1);
    }

    free(orientation1.adjacencyList);
    free(orientation1.reverseAdjacencyList);
    free(orientation2.adjacencyList);
    free(orientation2.reverseAdjacencyList);
}


int main(int argc, char ** argv) {
    struct options options = {.bruteForceFlag = false, .complementFlag = false,
     .exhaustiveCheckFlag = true, .doublecheckFlag=false,
     .oddCyclesHeuristicFlag = true, .verboseFlag = false, .printFlag = false, 
     .singleGraphFlag = false, .modulo = 1, .remainder = 0, 
     .sizeOfArray = 100000};
    struct counters numberOf = {0};
    int opt;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = 
        {   
            {"only-heuristic", no_argument, NULL, '2'},
            {"brute-force", no_argument, NULL, 'b'},
            {"complement", no_argument, NULL, 'c'},
            {"double-check", no_argument, NULL, 'd'},
            {"only-exact", no_argument, NULL, 'e'},
            {"help", no_argument, NULL, 'h'},
            {"print-orientation", no_argument, NULL, 'p'},
            {"single-graph-parallel", no_argument, NULL, 's'},
            {"verbose", no_argument, NULL, 'v'}
        };

        opt = getopt_long(argc, argv, "2bcdehpsv", long_options, &option_index);
        if (opt == -1) break;
        switch(opt) {
            case '2':
                options.exhaustiveCheckFlag = false;
                fprintf(stderr,
                 "Warning: fn can still be 2 even if output says >= 3.\n");
                fprintf(stderr, "Only using heuristic method.\n");
                break;
            case 'b':
                fprintf(stderr,
                 "Using brute force method where an exact method is used.\n");
                options.bruteForceFlag = true;
                break;
            case 'c':
                options.complementFlag = true;
                break;
            case 'd':
                options.doublecheckFlag = true;
                break;
            case 'e':
                fprintf(stderr, "Only using exact method.\n");
                options.oddCyclesHeuristicFlag = false;
                break;
            case 'h':
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case 'p':
                options.printFlag = true;
                options.verboseFlag = true;
                break;
            case 's':
                options.singleGraphFlag = true;
                break;
            case 'v':
                options.verboseFlag = true;
                break;
            case '?':
                fprintf(stderr,"Error: Unknown option: %c\n", optopt);
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr,
                 "Use ./findFrankNumber --help for more detailed instructions.\n");
                return 1;
        }
    }

    //  Check for other non-option arguments.
    bool haveModResPair = false;
    while (optind < argc) {
        bool pairIsInvalid = false;
        char* endptr;
        options.remainder = strtol(argv[optind], &endptr, 10);
        if( !endptr || *endptr != '/' || *(endptr+1) == '\0') {
            pairIsInvalid = true;
        }
        options.modulo = strtol(endptr+1, &endptr, 10);
        if( !endptr || *endptr != '\0') {
            pairIsInvalid = true;
        }
        if(options.modulo <= options.remainder) {
            pairIsInvalid = true;
        }
        if(haveModResPair) {
            fprintf(stderr,
             "Error: You can only add one res/mod pair as an argument.\n");
            fprintf(stderr, "%s\n", USAGE);
            fprintf(stderr,
             "Use ./findFrankNumber --help for more detailed instructions.\n");
            return 1;
        }
        if(pairIsInvalid) {
            fprintf(stderr,
                 "Error: Invalid res/mod pair: '%s'.\n", argv[optind]);
            fprintf(stderr, "%s\n", USAGE);
            fprintf(stderr,
             "Use ./findFrankNumber --help for more detailed instructions.\n");
            return 1;
        }
        fprintf(stderr, "Class=%d/%d.\n", options.remainder, options.modulo);
        haveModResPair = true;
        optind++;
    }
    if(options.oddCyclesHeuristicFlag) {
        fprintf(stderr,
         "Warning: this only works for cyclically 4-edge-connected graphs!\n");
    }
    if(options.printFlag && options.bruteForceFlag) {
        options.printFlag = false;
        fprintf(stderr,
         "Warning: no orientations will be printed for the brute force method.\n");
    }

    fprintf(stderr, "%s\n", 
     "Assuming graphs to be cubic and 3-edge-connected.");
    
    unsigned long long int totalGraphs = 0;
    unsigned long long int counter = 0;
    unsigned long long int skippedGraphs = 0;
    unsigned long long int passedGraphs = 0;
    clock_t start = clock();

    //  Start looping over lines of stdin.
    char * graphString = NULL;
    size_t size;
    while(getline(&graphString, &size, stdin) != -1) {
        totalGraphs++;
        numberOf.generatedOrientations = 0;
        numberOf.orientationsGivingSubset = 0;
        numberOf.complementaryBitsets = 0;
        numberOf.emptyBitsetsStored = 0;

        if(options.singleGraphFlag && totalGraphs >= 2) {
            fprintf(stderr, "Warning: do not input two graphs with -s.\n");
            totalGraphs--;
            break;
        }

        //  Skip graphs not belonging to res/mod class if singleGraphFlag is not
        //  active.
        if(!options.singleGraphFlag && 
         (totalGraphs - 1) % options.modulo != options.remainder) {
            continue;
        }

        int numberOfVertices = getNumberOfVertices(graphString);
        if(numberOfVertices == -1 || numberOfVertices > MAXVERTICES) {
            if(options.verboseFlag){
                fprintf(stderr, "Skipping invalid graph!\n");
            }
            skippedGraphs++;
            continue;
        }

        //  MAXVERTICES also indicates the largest size of a bitset, since we
        //  store edges in a bitset, the number of edges in a cubic graph
        //  (3*n/2) may not exceed MAXVERTICES.
        if(numberOfVertices*3/2 > MAXVERTICES) {
            if(options.verboseFlag){
                fprintf(stderr, "Skipping invalid graph! Too many edges.\n");
            }
            skippedGraphs++;
            continue;
        }
        bitset adjacencyList[numberOfVertices];
        if(loadGraph(graphString, numberOfVertices, adjacencyList) == -1) {
            if(options.verboseFlag){
                fprintf(stderr, "Skipping invalid graph!\n");
            }
            skippedGraphs++;
            continue;
        }
        counter++;


        if(options.verboseFlag) {
            fprintf(stderr, "Looking at:\n%s", graphString);
        }

        if(options.printFlag) {
            fprintf(stderr, "Labelling of graph:\n");
            printGraph(adjacencyList, numberOfVertices);
        }

        int frankNumber = 0;
        if(options.oddCyclesHeuristicFlag) {
            int F[numberOfVertices];
            if(hasSufficientCondition(adjacencyList, numberOfVertices, &options,
             &numberOf, complement(EMPTY, numberOfVertices), F)) {
                numberOf.graphsSatisfyingOddnessCondition++;
                frankNumber = 2;
            }
            else {
                if(options.verboseFlag) {
                    fprintf(stderr, 
                     "\tHeuristic failed. %soing exhaustive check.\n",
                     options.exhaustiveCheckFlag ? "D" : "Not d");
                }
                numberOf.graphsNotSatisfyingOddnessCondition++;
            }
        }
        if(options.exhaustiveCheckFlag && frankNumber == 0) {
            frankNumber = findFrankNumber(adjacencyList, numberOfVertices, 
                &options, &numberOf);
            if(options.verboseFlag) {
                fprintf(stderr,
                 "\tStrongly connected orientations generated: %llu\n",
                 numberOf.generatedOrientations);
                if(options.bruteForceFlag) {
                    fprintf(stderr, "\tOrientations giving subsets: %llu\n",
                     numberOf.orientationsGivingSubset);
                    fprintf(stderr, "\tOrientations giving supersets: %llu\n",
                     numberOf.orientationsGivingSuperset);
                    fprintf(stderr, "\tNumberOfComplementaryBitsets: %llu\n",
                     numberOf.complementaryBitsets);
                }
            }
        }
        if(frankNumber == 0) {
            if(options.verboseFlag) {
                fprintf(stderr, "\tFrankNumber >= 3.\n\n");
                fprintf(stderr, "------------------------------------\n\n");
            }
            if(!options.complementFlag) {
                passedGraphs++;
                printf("%s", graphString);
            }
        }
        if(frankNumber == 2) {
            if(options.verboseFlag) {
                fprintf(stderr, "\tFrankNumber = 2.\n\n");
                fprintf(stderr, "------------------------------------\n\n");
            }
            if(options.complementFlag) {
                passedGraphs++;
                printf("%s", graphString);
            }
        }
        if(numberOf.mostGeneratedOrientations < numberOf.generatedOrientations) {
            numberOf.mostGeneratedOrientations = numberOf.generatedOrientations;
        }
        if(numberOf.mostStoredBitsets < numberOf.storedBitsets) {
            numberOf.mostStoredBitsets = numberOf.storedBitsets;
        }

    }
    free(graphString);
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    if(options.bruteForceFlag) {
        fprintf(stderr, 
         "Largest size of bitset array is %llu elements (%.2f GB)\n",
          numberOf.mostStoredBitsets, numberOf.mostStoredBitsets*8/1000000000.0);
    }
    fprintf(stderr,"\rChecked %lld graphs in %f seconds: %llu %s.\n",
     counter, time_spent, passedGraphs, options.complementFlag ? 
     (options.exhaustiveCheckFlag ? 
     "have fn = 2": 
     "passed sufficient condition for fn 2") : 
     (options.exhaustiveCheckFlag ? "have fn > 2" : 
     "did not pass sufficient condition for fn 2"));
    if(skippedGraphs > 0) {
        fprintf(stderr, "Warning: %lld graphs were skipped.\n", skippedGraphs);
    }
    if(options.oddCyclesHeuristicFlag) {
        fprintf(stderr, 
         "%llu satisfied at least one of the sufficient conditions. %llu did not.\n", 
         numberOf.graphsSatisfyingOddnessCondition,
         numberOf.graphsNotSatisfyingOddnessCondition);
        fprintf(stderr, "%llu satisfied first and %llu satisfied second\n",
         numberOf.graphsSatisfyingFirstOddness,
         numberOf.graphsSatisfyingSecondOddness);
    }

    return 0;
}