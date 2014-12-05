// Class: Graph Code

#include "Graph.h"

struct Graph::node {
	node() { }
	node(int _ID, Sequence* _seq) : ID(_ID), seq(_seq) { }
	int ID;
	Sequence* seq;
	edge* child[2];
};

struct Graph::edge {
	edge() { }
	edge(int _to, double _weight) : to(_to), weight(_weight) { }
	int to;
	double weight;
};
	
Graph::Graph() { }

Graph::~Graph() { }

void Graph::addLeafNode(Sequence* seq) {
	// creates a node for an individual sequence
	G.push_back(new node(static_cast<int>(G.size()), seq));
}

void Graph::addEdge(int from, int to, double weight) {
	// needed? Probably not

}

void Graph::addParentNode(node* child1, node* child2) {
	// create new node with edges to the children
	int index = G.size();
	double weight = 0;
	// Align sets of sequences here?
	/*
	*	We can try this...
	*	We assume each set of sequences is locked in step. Meaning each index of the sequences maintains an array of characters for each sequence in that alignment.
	*	When we go to align the two sets, we look at each index and compare the two sets of char at each site. We then introduce gaps into a set of aligned sequences. Scoring would be based on the total number of chars. So, an index where all are the same could recieve a higher score than an index where many are different
	*	Maybe use the sum of all pairs at an index to score that position
	*	EX:
			A  A  T  G  C  A
			A  A  T  G  G  A
			A  A  T  G  T  A

			A  A  T  G  T  A
			A  A  T  G  T  A

	score:	12 12 12 12 0  12

	*	So, we could use a modified NW Align algorithm, where the score at a particular pair of indices is the sum of all the pairs of char between the two alignments at that index.
	*/
	G.push_back(new node(index, NULL));
	G[index]->child[0] = new edge(child1->ID, weight);
	G[index]->child[1] = new edge(child2->ID, weight);


}
