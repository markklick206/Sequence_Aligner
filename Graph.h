// Class: Graph

#include <iostream>
#include <vector>

#include "Sequence.h"

class Graph {
public:
	struct node;
	struct edge;
	
	Graph();
	~Graph();
	
	void addLeafNode(int edgeval);
	void addEdge(int to, int from);
	void joinNodes(node* node1, node* node2, node* parent);
	
	
private:
	node* head;
	vector<node*> G;
};
