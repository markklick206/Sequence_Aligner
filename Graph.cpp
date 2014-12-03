// Class: Graph

#include "Graph.h"

struct node {
	vector<edge*> edges;
	int ID;
	Sequence seq;
};

struct edge {
	node* to;
	double weight;
};
	
Graph::Graph() { head = 0; }

Graph::~Graph() { }

void Graph::addLeafNode(int edgeVal, int from) {
	//creates a node for each individual sequence
	
}

void Graph::addEdge(int to, int from) {
	//needed?

}

void Graph::addParentNode(node* child1, node* child2) {
	//create new node with edges to the children
}
