#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>

#include <GfaGraph.h>
#include "gfa2graph.h"
#include "graph.h"
#include "io.h"

GfaGraph load_gfa(const std::string &gfa_filename) {
	GfaGraph gfa = GfaGraph::LoadFromFile(gfa_filename);
	LOG_INFO << "GFA loaded with " << gfa.nodes.size() << " nodes and " << gfa.edges.size() << " edges";
	for (const auto edge: gfa.edges) {
		assert(edge.first.end);
		for (const auto n: edge.second) {
			assert(n.end);
		}
	}
	assert(gfa.edgeOverlap <= 0);
	return gfa;
}

void gfa2graph(GfaGraph &gfa, astarix::graph_t *G) {
	LOG_INFO << "GFA to Internal graph";

	std::unordered_map<int, int> node2idx;

	for (const auto &node: gfa.nodes) {
		node2idx[node.first] = G->add_node();
		assert(!node.second.empty());
	}

	for (const auto &node: gfa.nodes) {
		int curr = node2idx[node.first];
		for (int j=0; j<node.second.size()-1; j++) {
			int next = G->add_node();
			G->add_edge(curr, next, node.second[j], astarix::ORIG);
			curr = next;
		}
		auto last_letter = node.second.back();
		for (const auto &neigh: gfa.edges[NodePos(node.first, '+')]) {
			G->add_edge(curr, node2idx[neigh.id], last_letter, astarix::ORIG);
		}
	}

	G->orig_nodes = G->V.size();
	G->orig_edges = G->E.size();
}
