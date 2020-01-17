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
	int cnt;
	TrieNode* children[4];

	TrieNode(int _node)
		: node(_node), cnt(0) {
		children[0] = children[1] = children[2] = children[3] = nullptr;
	}

	TrieNode* get_node(char label) {
		++cnt;
		auto &p = children[nucl2num(label)];
		if (p == nullptr)
			p = new TrieNode(-1);
		return p;
	}

	void del_node() {
		for (int i=0; i<4; i++) {
			if (children[i] != nullptr) {
				children[i]->del_node();
				delete children[i];
			}
		}
	}
};

void dfs_construct_trie(const graph_t &G, int v, int rem_depth, TrieNode *tree_v) {
	for (int idx=G.V[v]; idx!=-1; idx=G.E[idx].next) {
		const edge_t &e = G.E[idx];
		assert(e.type == ORIG);
		if (rem_depth > 0) {
			TrieNode *next_tree_v = tree_v->get_node(e.label);
			dfs_construct_trie(G, e.to, rem_depth-1, next_tree_v);
		}
	}
}

void dfs_trie_to_graph(const graph_t &G, int v, int rem_depth, TrieNode *tree_v, EdgeList *new_edges, int *curr_node) {
	if (tree_v->node == -1) {
		tree_v->node = *curr_node;
		++(*curr_node);
	}

	for (int idx=G.V[v]; idx!=-1; idx=G.E[idx].next) {
		const edge_t &e = G.E[idx];
		if (rem_depth == 0 || tree_v->cnt == 1) {
			// Connect to reference genome.
			new_edges->push_back(make_pair(tree_v->node, make_pair(e.to, e.label)));
		} else {
			// Connect to deeper trie node.
			TrieNode *next_tree_v = tree_v->get_node(e.label);
			bool new_edge = (next_tree_v->node == -1);
			dfs_trie_to_graph(G, e.to, rem_depth-1, next_tree_v, new_edges, curr_node);
			if (new_edge)
				new_edges->push_back(make_pair(tree_v->node, make_pair(next_tree_v->node, e.label)));
		}
	}
}

void add_tree(graph_t *G, int tree_depth) {
	TrieNode tree_root(0);
	EdgeList new_edges;
	int curr_node=G->V.size();
	G->trie_first_node = curr_node;

	try {
		// Construct Trie.
		for (int i=1; i<G->V.size(); i++)
			dfs_construct_trie(*G, i, tree_depth, &tree_root);
		// Extrect edges to be added to the graph.
		for (int i=1; i<G->V.size(); i++)
			dfs_trie_to_graph(*G, i, tree_depth, &tree_root, &new_edges, &curr_node);
	} catch (std::bad_alloc& ba) {
		std::cerr << "new_edges.size(): " << new_edges.size() << '\n';
		std::cerr << "bad_alloc caught: " << ba.what() << '\n';
		throw;
    }

	G->trie_nodes = curr_node - G->V.size();
	G->trie_edges = new_edges.size();

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

	tree_root.del_node();
}
