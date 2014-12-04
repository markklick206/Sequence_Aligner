// Class: Graph Header

#include <iostream>
#include <vector>

#include "Sequence.h"

class Graph {
public:
	struct node;
	struct edge;
	
	Graph();
	~Graph();
	
	void addLeafNode(Sequence* seq);
	void addEdge(int from, int to, double weight);
	void addParentNode(node* child1, node* child2);
	
private:
	std::vector<node*> G;
};
