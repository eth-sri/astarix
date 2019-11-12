#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>

#include "graph.h"
#include "io.h"

using namespace std;

int main(int argc, char **argv) {
	astarix::init_logger();

	if (argc != 2) {
		cerr << "Fasta file name expected" << endl;
		return 1;
	}

	// Read fasta.
	vector<astarix::seq_t> fastas = astarix::read_fasta(argv[1]);
	assert(fastas.size() == 1);
	astarix::seq_t fasta = *fastas.begin();

	// Build graph.
	astarix::graph_t G;
	int source=G.add_node();
	int sink=G.add_node();
	G.add_seq(1.0, source, fasta.s, sink);
	G.write(stdout);

	return 0;
}
