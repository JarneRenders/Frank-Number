/**
 *  This header file contains functions that have to do with checking
 *  hamiltonicity properties of graphs (K1-hamiltonicity, K2-hamiltonicity).
 * */

#ifndef HAM_METHODS 
#define HAM_METHODS

#include "bitset.h"

/**
 *  Returns a boolean indicating whether or not the specified path can be
 *  extended to a hamiltonian cycle in the specified graph. The path is
 *  represented by its first and last element, length and the absence of its
 *  elements in remainingVertices.
 *  
 *  @param  adjacencyList   An array of bitsets representing the adjacency
 *   list of the graph. We check whether the path can be extended to a
 *   hamiltonian cycle in this graph.
 *  @param  remainingVertices   A bitset containing all vertices of the graph
 *   we still consider adding to the path. One can also remove vertices from
 *   this bitset on initialization to ban vertices from ever appearing in the
 *   path.
 *  @param  lastElemOfPath  An int representing the last vertex in the path.
 *  @param  firstElemOfPath An int representing the first vertex in the path.
 *  @param  numberOfVertices    The number of vertices in the graph.
 *  @param  pathLength  The length of the current path.
 * 
 *  @return Boolean representing whether the given path can be extended to a
 *   hamiltonian cycle.
 * */
bool canBeHamiltonian(bitset adjacencyList[], bitset remainingVertices, int
lastElemOfPath, int firstElemOfPath, int numberOfVertices, int pathLength);


/**
 * Similar to the canBeHamiltonian, but specifically for counting and printing
 * cycles/paths. Has slightly worse performance than canBeHamiltonian.
 * 
 *  @param  adjacencyList   An array of bitsets representing the adjacency
 *   list of the graph. We check whether the path can be extended to a
 *   hamiltonian cycle in this graph.
 *  @param  remainingVertices   A bitset containing all vertices of the graph
 *   we still consider adding to the path. One can also remove vertices from
 *   this bitset on initialization to ban vertices from ever appearing in the
 *   path.
 *  @param  pathList    A list of ints representing the current path.
 *  @param  lastElemOfPath  An int representing the last vertex in the path.
 *  @param  firstElemOfPath An int representing the first vertex in the path.
 *  @param  numberOfVertices    The number of vertices in the graph.
 *  @param  pathLength  The length of the current path.
 *  @param  numberOfhamiltonianCycles   Pointer to an int representing
 *   counting the number of hamiltonian cycles in the graph. 
 *  @param  allCyclesFlag   Boolean which if true will count the number of
 *   hamiltonian cycles containing the starting path.
 *  @param  verboseFlag Boolean which if true will print a hamiltonian cycle
 *   if one is found or all hamiltonian cycles if allCyclesFlag is true/
 * 
 *  @return Boolean representing whether the given path can be extended to a
 *   hamiltonian cycle.
 * 
 */
bool canBeHamiltonianPrintCycle(bitset adjacencyList[], bitset
remainingVertices, int pathList[], int lastElemOfPath, int firstElemOfPath,
int numberOfVertices, int pathLength, int* numberOfhamiltonianCycles, bool
allCyclesFlag, bool verboseFlag);

/**
 *  Returns a boolean indicating whether the subgraph of the given graph
 *  spanned by all vertices not in excludedVertices is hamiltonian or not. It
 *  does so by checking whether all paths of the form a, v, b, where v is a
 *  vertex not in excludedVertices of lowest degree, a, b are neighbours of v
 *  not in excludedVertices and a > b, can be extended to some hamiltonian
 *  cycle.
 * 
 *  @param  adjacencyList   An array of bitsets representing the adjacency
 *   list of the given graph.
 *  @param  numberOfVertices    The number of vertices in the graph.
 *  @param  excludedVertices    A bitset representing the vertices we ban from
 *   our graph. We check whether the subgraph spanned by all vertices not
 *   contained in exludedVertices is hamiltonian or not.
 *  @param  allCyclesFlag   Boolean which indicates whether we compute and
 *   print out the number of cycles in the (sub)graph.
 *  @param verboseFlag  Boolean indication whether we should print out a
 *   hamiltonian cycle if it is found or all hamiltonian cycles if
 *   allCyclesFlag is true.
 * 
 *  @return Returns true when the (sub)graph is hamiltonian and false
 *   otherwise.
 * */
bool isHamiltonian(bitset adjacencyList[], int numberOfVertices, bitset
excludedvertices, bool allCyclesFlag, bool verboseFlag);


/**
 * Returns a boolean indicating whether or not a given graphs has
 * minimumdegree degree.
 * 
 * @param   adjacencyList   An array of bitsets representing the adjacency
 *  list of the given graph.
 * @param   numberOfVertices    The number of vertices in the graph.
 * @param   degree  The degree we want to check.
 * 
 * @return Returns true when the graph has minimumdegree degree and false
 *  otherwise.
 */
bool hasMinimumDegree(bitset adjacencyList[], int numberOfVertices, int degree);

/**
 *  Returns a boolean indicating whether the graph is K1-hamiltonian, i.e.
 *  deleting any copy of K1 (a single vertex), yields a hamiltonian graph for
 *  every vertex.
 * 
 *  @param  adjacencyList   Array of bitsets representing the adjacency list
 *   of the original graph.
 *  @param  numberOfVertices The number of vertices in the original graph.
 *  @param  verboseFlag     Boolean, which if true indicates which vertex-deleted
 *   subgraphs are non-hamiltonian (i.e. which vertices are exceptional). If
 *   vertexToCheck is in the graph it prints out (if any exist) a hamiltonian
 *   cycle in G - vertexToCheck.
 *  @param  allCyclesFlag   Boolean, which if true and if vertexToCheck is in
 *   the graph counts or prints out all hamiltonian cycles in G -
 *   vertexToCheck depending on the value of verboseFlag.
 *  @param  vertexToCheck   If this represents a vertex in the graph a
 *   hamiltonian cycle (if any exist) of G - vertexToCheck will be printed if
 *   verboseFlag is true or all hamiltonian cycles will be printed if
 *   allCyclesFlag and verboseFlag are true or all hamiltonian cycles will be
 *   counted if allCyclesFlag is true, but verboseFlag is false.  
 * 
 *  @return True if the graph is K1-hamiltonian.
 * */
bool isK1Hamiltonian(bitset adjacencyList[], int numberOfVertices, bool
verboseFlag, bool allCyclesFlag, int vertexToCheck);

/**
 *  Returns a boolean indicating whether the graph is K2-hamiltonian, i.e.
 *  deleting any copy of K2 (a pair of adjacent vertices), yields a
 *  hamiltonian graph for all such copies.
 * 
 *  @param  adjacencyList   Array of bitsets representing the adjacency list
 *   of the original graph.
 *  @param  numberOfVertices The number of vertices in the original graph.
 *  @param  verboseFlag     Boolean, which if true indicates which
 *   adjacent-pair-deleted subgraphs are non-hamiltonian. If
 *   vertexPairToCheck is a pair of adjacent vertices in the graph, a
 *   hamiltonian cycle (if any exist) in G - vertexPairToCheck
 *   [0] - vertexPairToCheck[1] will be printed.
 *  @param  allCyclesFlag   Boolean, which if true and if vertexPairToCheck is
 *   in the graph counts or prints out all hamiltonian cycles in G -
 *   vertexPairToCheck[0] - vertexPairToCheck[1] depending on the value of
 *   verboseFlag.
 *  @param  vertexPairToCheck   List representing a pair of vertices. If these
 *   are adjacent vertices in the graph a hamiltonian cycle (if any exist) of
 *   G - vertexPairToCheck[0] - vertexPairToCheck[1] will be printed if
 *   verboseFlag is true or all hamiltonian cycles will be printed if
 *   allCyclesFlag and verboseFlag are true or all hamiltonian cycles will be
 *   counted if allCyclesFlag is true, but verboseFlag is false.  
 * 
 *  @return True if the graph is K2-hamiltonian.
 * */
bool isK2Hamiltonian(bitset adjacencyList[], int numberOfVertices, bool
verboseFlag, bool allCyclesFlag, int vertexPairToCheck[]);


/**
 * Returns an integer indicating whether or not the (sub)graph contains a
 * hamiltonian path with specified endpoints.
 * 
 * @param  adjacencyList   An array of bitsets representing the adjacency list
 *  of the given graph.
 * @param  numberOfVertices    The number of vertices in the graph.
 * @param  excludedVertices    A bitset representing the vertices we ban from
 *  our graph. We check whether the subgraph spanned by all vertices not
 *  contained in exludedVertices contains a hamiltonian path or not.
 * @param   start   One of the endpoints between which we want to find a
 *  hamiltonian path.
 * @param   end     The other endpoint between which we want to find a
 *  hamiltonian path.
 * @param   verboseFlag     If this boolean is true, we print out a
 *  hamiltonian path between start and end, if one exists.
 * @param   allCyclesFlag   If this boolean is true, we count how many
 *  hamiltonian paths there are between start and end. If verboseFlag is
 *  true, we print all of these.
 * 
 * @return Non-zero if there exists a hamiltonian path between start and end
 *  in the (sub)graph and zero otherwise. If allCyclesFlag is present, the
 *  return value is the number of hamiltonian paths between the endpoints
 * */
int containsHamiltonianPathWithEnds(bitset adjacencyList[], int
numberOfVertices, bitset excludedNodes, int start, int end, bool verboseFlag,
bool allCyclesFlag);

/**
 * Returns a boolean indicating whether or not the (sub)graph contains two
 * disjoint paths with specified endpoints and containing specified vertices,
 * which together span the (sub)graph.
 * 
 *  @param  adjacencyList   An array of bitsets representing the adjacency
 *   list of the given graph.
 *  @param  numberOfVertices    The number of vertices in the graph.
 *  @param  excludedVertices    A bitset representing the vertices we ban from
 *   our graph. We check whether the subgraph spanned by all vertices not
 *   contained in exludedVertices contains two disjoint spanning paths.
 *  @param  startOfPath1   One of the endpoints of the first path.
 *  @param  endOfPath1  The other endpoint of the first path.
 *  @param  verticesContainedByPath1    A bitset of vertices we require to be
 *   in the first path.
 *  @param  startOfPath2    One of the endpoints of the second path.
 *  @param  endOfPath2  The other endpoint of the second path.
 *  @param  verticesContainedByPath2    A bitset of vertices we require to be
 *   in the second path.
 *  @param  allCyclesFlag   If this boolean is true, we count how many pairs
 *   of disjoint paths starting with the specified endpoints and containing
 *   the specified vertices span the (sub)graph. If verboseFlag is true all
 *   these paths get printed.
 *  @param verboseFlag  If this boolean is true, we print a pair of disjoint
 *   paths with specified endpoints and containing the specified vertices
 *   that span the graph (if such a pair exists).
 * */
bool containsDisjointSpanningPathsWithEnds(bitset adjacencyList[], int
numberOfVertices, bitset excludedVertices, int startOfPath1, int endOfPath1,
bitset verticesContainedByPath1, int startOfPath2, int endOfPath2, bitset
verticesContainedByPath2, bool allCyclesFlag, bool verboseFlag);


/**
 * Returns a boolean indicating whether or not the (sub)graph is traceable,
 * i.e., contains some hamiltonian path.
 *  @param  adjacencyList   An array of bitsets representing the adjacency
 *   list of the given graph.
 *  @param  numberOfVertices    The number of vertices in the graph.
 *  @param  excludedVertices    A bitset representing the vertices we ban from
 *   our graph. We check whether the subgraph spanned by all vertices not
 *   contained in exludedVertices contains two disjoint spanning paths.
 *  @param  allCyclesFlag   If this boolean is true, we count how many
 *   hamiltonian paths there are in the graph. If verboseflag is also true
 *   all these paths get printed and sorted by their endpoints.
 *  @param  verboseFlag If this boolean is true, we print a hamiltonian path.
 * 
 *  @return True if (sub)graph is traceable, false otherwise.
 * */
bool isTraceable(bitset adjacencyList[], int numberOfVertices, bitset
excludedVertices, bool allCyclesFlag, bool verboseFlag);

/**
 *  Returns a boolean indicating whether the graph is K1-traceable, i.e.
 *  deleting any copy of K1 (a single vertex), yields a traceable graph for
 *  every vertex.
 * 
 *  @param  adjacencyList   Array of bitsets representing the adjacency list
 *   of the original graph.
 *  @param  numberOfVertices The number of vertices in the original graph.
 *  @param  verboseFlag     Boolean, which if true indicates which vertex-deleted
 *   subgraphs are non-traceable (i.e. which vertices are exceptional). If
 *   vertexToCheck is in the graph it prints out (if any exist) a hamiltonian
 *   path in G - vertexToCheck.
 *  @param  allCyclesFlag   Boolean, which if true and if vertexToCheck is in
 *   the graph counts or prints out all hamiltonian paths in G -
 *   vertexToCheck depending on the value of verboseFlag.
 *  @param  vertexToCheck   If this represents a vertex in the graph a
 *   hamiltonian path (if any exist) of G - vertexToCheck will be printed if
 *   verboseFlag is true or all hamiltonian paths will be printed if
 *   allCyclesFlag and verboseFlag are true or all hamiltonian cycles will be
 *   counted if allCyclesFlag is true, but verboseFlag is false.  
 * 
 *  @return True if the graph is K1-traceable, false otherwise.
 * */
bool isK1Traceable(bitset adjacencyList[], int numberOfVertices, bool
allCyclesFlag, bool verboseFlag, int vertexToCheck);

#endif