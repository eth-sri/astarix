#include "graph.h"
#include "io.h"

using namespace std;
using namespace astarix;

graph_t G;

bool get_next_edge(node_t u, edge_t *e) {
	for (auto it=G.begin_orig_edges(u); it!=G.end_orig_edges(); ++it) {
		//cerr << e << endl;
		*e = *it;
		return true;
	}
	return false;
}

int main(int argc, char **argv) {
    ios_base::sync_with_stdio(false);

	if (argc != 2) {
		cerr << "One parameter is expected: a .gfa filename" << endl;
		return 1;
	}

	string graph_file = argv[1];
	GfaGraph gfa = load_gfa(graph_file);
	gfa2graph(gfa, &G);

	node_t u;  // Assuming 1 is "the beginning" of the genome
	for (u=1; u<G.nodes(); u++)
		if (!G.hasIncomingEdges(u))
			break;
	for (node_t v=u+1; v<G.nodes(); v++)
		if (!G.hasIncomingEdges(v)) {
			cerr << "Warning: There is at least one more path with a different start." << endl;
			break;
		}
			
	if (u == G.nodes())
		exit(1);
	// u does not have incoming edges

	// graph
	string L;
	edge_t e;
	while (get_next_edge(u, &e)) {
		L += e.label;
		u = e.to;
	}
	cerr << "gfa nodes = " << G.nodes() << endl;
	cerr << "fa length = " << L.size() << endl;
	
	// fasta
	cout << ">Created from " << graph_file << " using gfa2fa.cpp" << endl;
	cout << L << endl;

	return 0;
}
