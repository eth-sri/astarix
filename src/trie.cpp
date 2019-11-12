#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <string>

#include "graph.h"
#include "io.h"

using namespace std;
using namespace astarix;

typedef vector<pair<int, pair<int, char>>> EdgeList;

struct TrieNode {
	int node;
	map<char, TrieNode*> children;

	TrieNode(int _node)
		: node(_node) {}

	TrieNode* add_node(const graph_t &G, char label, EdgeList *new_edges, int *curr_node) {
		if (children.find(label) != children.end()) {
			return children[label];
		} else {
			int v = *curr_node;
			++(*curr_node);
			TrieNode* tree_v = new TrieNode(v);
			children[label] = tree_v;
			new_edges->push_back(make_pair(node, make_pair(v, label)));
			return tree_v;
		}
	}
};

void dfs(const graph_t &G, int v, int rem_depth, TrieNode *tree_v, EdgeList *new_edges, int *curr_node) {
	for (int idx=G.V[v]; idx!=-1; idx=G.E[idx].next) {
		const edge_t &e = G.E[idx];
		assert(e.type == ORIG);
		if (rem_depth > 0) {
			TrieNode *next_tree_v = tree_v->add_node(G, e.label, new_edges, curr_node);
			dfs(G, e.to, rem_depth-1, next_tree_v, new_edges, curr_node);
		} else {
			new_edges->push_back(make_pair(tree_v->node, make_pair(e.to, e.label)));
		}
	}
}

void add_tree(graph_t *G, int tree_depth) {
	TrieNode tree_root(0);
	EdgeList new_edges;
	int curr_node=G->V.size();

	try {
		// Construct Trie.
		for (int i=1; i<G->V.size(); i++)
			dfs(*G, i, tree_depth, &tree_root, &new_edges, &curr_node);
	} catch (std::bad_alloc& ba) {
		std::cerr << "new_edges.size(): " << new_edges.size() << '\n';
		std::cerr << "bad_alloc caught: " << ba.what() << '\n';
		throw;
    }

	try {
		// Add the tree edges to astarix
		while (G->V.size() < curr_node)
			G->add_node();
		for (const auto &e: new_edges) {
			G->add_edge(e.first, e.second.first, e.second.second, JUMP);
		}
	} catch (std::bad_alloc& ba) {
		std::cerr << "G->V.size(): " << G->V.size() << '\n';
		std::cerr << "bad_alloc caught: " << ba.what() << '\n';
		throw;
    }
}
