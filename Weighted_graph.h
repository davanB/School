/*****************************************
 * Code by Davan Basran
 * Prepared as a requirement for school
 *****************************************/

#ifndef WEIGHTED_GRAPH_H
#define WEIGHTED_GRAPH_H

#ifndef nullptr
#define nullptr 0
#endif

#include <iostream>
#include <limits>
#include "Exception.h"
#include "Disjoint_sets.h"
#include <iterator>
#include <set>

/*
class Weighted_graph
used to hold a weighted graph to implement Kruskals Algorithm

Private members:
	static const double INF;
	int size; // num of verticies
	double **graph_array; // array of pointers act as the 'columns' of the adjacency matrix
	int *indeg; // array keeps track of number of edges verticie has
	int edges_counted; //keeps track of number of edges in graph

Public Accessors:
	int degree( int ) const;
	returns the number of verticies that are adjacent to the passed vertex; or the number of edges it has

	int edge_count() const;
	returns the total number of edges in the graph
*/
class Weighted_graph {
	private:
		static const double INF;
		int size;
		double **graph_array;
		int *indeg; 
		int edges_counted; 

		Weighted_graph( Weighted_graph const & );
		Weighted_graph &operator=( Weighted_graph );

		// your choice

	public:
		Weighted_graph( int = 10 );
		~Weighted_graph();

		int degree( int ) const;
		int edge_count() const;
		std::pair<double, int> minimum_spanning_tree() const;

		bool insert_edge( int, int, double );
		bool erase_edge( int, int );
		void clear_edges(); 

	// Friends

	friend std::ostream &operator<<( std::ostream &, Weighted_graph const & );
};

// used to neatly packages a pair of verticies that form an edge
typedef std::pair<int, int> vertex;

// used to store an edge and all relavent informaiton about it, a pair of verticies and its weight
struct edge 
{
	vertex v;
	double weight;
};

// used by the std::set to compare the weights of the edges it stores and sort them accoordingly
struct compare
{
	bool operator()(const edge &e, const edge &e2) {
		return e.weight < e2.weight;
	}
};

// initialize a std::set to store and sort edges based on weight and an iterator to iterate through the set and access the edges
std::multiset<edge, compare>  edges;
std::multiset<edge, compare>::iterator it;

// global cont INF
const double Weighted_graph::INF = std::numeric_limits<double>::infinity();

/*
Weighted_graph::Weighted_graph(int n)
Constuctor for class Weighted_graph

creates an upper triangle adjacency matrix to store the weights of edges between 2 verticies
All other entries are given INF weight to signify they are not connected

Uses pointer arithmetic to make the pointer of pointers act as the columns of the matrix thus can access entires like for a matrix
all entires are actually stored in an array

An array indeg is used to keep track of the in degree of each vertex and is intialized and given intial 
value of 0 while we initalize all entires to have INF weight
*/
Weighted_graph::Weighted_graph(int n) :
size(std::max(1, n)), graph_array(new double* [size]), 
indeg(new int [size]), edges_counted(0) {
	// create the upper triangle matrix
	graph_array[0] = nullptr;
	graph_array[1] = new double[size*(size-1)/2];
	for (int i = 2; i < size; ++i)
	{
		graph_array[i] = graph_array[i - 1] + i -1;
	}
	for (int j = 0; j < size; j++)
	{
		for (int k = 0; k < j; k++)
		{
			graph_array[j][k] = INF;
		}
		indeg[j] = 0; // initialize in deg of all verticies
	}
}

/*
Weighted_graph::~Weighted_graph()
Destructor for weighted_graph
deallocated memory for all arrays but also clears the set of all data relavent to the deleted graph
*/
Weighted_graph::~Weighted_graph() {
	delete[] graph_array[1];
	delete[] graph_array;
	delete[] indeg;
	edges.clear();
}

// see comments for class for details
int Weighted_graph::degree(int i) const {
	if (i < 0 || i >(size - 1)) throw illegal_argument();
	return indeg[i];
}

//see class for details
int Weighted_graph::edge_count() const {
	return edges_counted;
}

/*
bool Weighted_graph::insert_edge( int i, int j, double d )
inserts edges into the graph

Param1:
	int i, which is the first vertex or the row in the matrix
Param 2:
	int j, which is the second vertex or the column in the matrix;
Param 3:
	double d; this is the weight of the edge which is to be inserted

There is no exception handling but there is error checking which throw the appropirate exceptions

Cases:
	if i == j return false as we cant insert an edge between the same point creating a loop
	if j < i then use i as the row and j as the column to reduce redundancy of having to update both
		If there was already an edge inserted there erase the entry in the set and lower the indeg count as we will be adding to it later again
		else insert a new entry in the set and update the matrix along with the indeg array and total number of edges
	Else do everything menioned above but simply make j the row and i the column
*/
bool Weighted_graph::insert_edge( int i, int j, double d ) {
	if (i < 0 || j < 0 || d < 0) throw illegal_argument();
	else if (i > size-1 || j > size-1) throw illegal_argument();
	else if (i == j) return false;
	else {
		if (j < i)  {
			edge e;
			vertex v;
			v.first = i;
			v.second = j;
			e.v = v;
			if (graph_array[i][j] != INF) { //we must delete whats in set before inserting again
				e.weight = graph_array[i][j];
				indeg[i]--;                  // to make sure double counting does not occur
				indeg[j]--;
				edges_counted--;
				edges.erase(e);
			}
			e.weight = d;
			graph_array[i][j] = d;
			edges.insert(e);
		}
		else {
			edge e;
			vertex v;
			v.first = i;
			v.second = j;
			e.v = v;
			if (graph_array[j][i] != INF) {
				e.weight = graph_array[j][i];
				indeg[i]--; // same steps as above
				indeg[j]--;
				edges_counted--;
				edges.erase(e);
			}
			e.weight = d;
			graph_array[j][i] = d;
			edges.insert(e);
		}
		indeg[i]++;
		indeg[j]++;
		edges_counted++;
		return true;
	}
}

/*
bool Weighted_graph::erase_edge(int i, int j)
used to remove an edge from the graph

Param 1:
	int i, which is one of the 2 vertex betweent the edge to be removed
Param 2:
	int j, which is the second vertex adacent to i which forms the edge to be removed

There is no exception handling but there is error checking to make sure paramaters are correct

Cases:
	If there is no edge where we expect there to be one then return false
	Else remove the corresponding entry from the set and from the graph matrix and update all relavent fields
*/
bool Weighted_graph::erase_edge(int i, int j) {
	if (i < 0 || j < 0) throw illegal_argument();
	else if (i > size-1 || j > size-1) throw illegal_argument();
	else if (i == j) return true;
	else {
		if (j < i)
			if (graph_array[i][j] == INF) return false; // edge does not exist
			else {
				vertex v;
				v.first = i; v.second = j;
				edge e;
				e.v = v;
				e.weight = graph_array[i][j];
				edges.erase(e);
				graph_array[i][j] = INF;
				
			}
		else
			if (graph_array[j][i] == INF) return false; //edge does not exist
			else {
				vertex v;
				v.first = j; v.second = i;
				edge e;
				e.v = v;
				e.weight = graph_array[j][i];
				edges.erase(e);
				graph_array[j][i] = INF;
			}
		edges_counted--;
		indeg[i]--;
		indeg[j]--;
		return true;
	}
}

/*
void Weighted_graph::clear_edges()
manually iterate through the graph and remove all edges and update all relavent fields
clear the set which contains all the graphs data 
*/
void Weighted_graph::clear_edges() {
	for (int j = 0; j < size; j++)
	{
		for (int k = 0; k < j; k++)
		{
			if (j == k) graph_array[j][k] = 0;
			else {
				graph_array[j][k] = INF;
			}
		}
		indeg[j] = 0;
	}
	edges.clear();
	edges_counted = 0;
}

/*
std::pair<double, int> Weighted_graph::minimum_spanning_tree() const
implements Kruskals Algorithm for the weighted graph constucted

There is no error checking or exception handling

We being by intializing a disjoint set to keep track of the verticies currently in the mst
Then starting with the first element in the set (the min entry) while we have no made the mst or there are no more edges left
	Check if the edges vertecies are in the mst or not (if they are in the same set)
		If not then take the union of them (combine the 2 sets) and add the weight of the edge to the total mst weight
	move on to the next element in the set

Finish by returning an std::pair of the weight and total number of edges we had to check to form the mst
*/
std::pair<double, int> Weighted_graph::minimum_spanning_tree() const {

	Data_structures::Disjoint_sets edge_sets(size); // make a disjoint set to keep track of verticies in min spanning tree
	int num_of_edges_checked = 0;
	int num = 0;
	double weight_of_tree = 0.0;
	it = edges.begin(); // start at the first element in the set
	while (num != (size - 1) && it != edges.end()) //until spanning tree made or until no more edges to check
	{
		num_of_edges_checked++;
		if (edge_sets.find(it->v.first) != edge_sets.find(it->v.second)) {
			edge_sets.set_union(it->v.first, it->v.second);
			num++;
			weight_of_tree += it->weight;
		}
		it++;
	}

	return std::pair<double, int>( weight_of_tree, num_of_edges_checked );
}

std::ostream &operator<<( std::ostream &out, Weighted_graph const &graph ) {
	for (int j = 0; j < graph.size; j++)
	{
		for (int k = 0; k < j; k++)
		{
			std::cout << '[' << graph.graph_array[j][k] << ']';
		}
		std::cout << std::endl;
	}
	for (int l = 0; l < graph.size; l++)
	{
		std::cout << graph.indeg[l];
	}

	return out;
}

#endif
