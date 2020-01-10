#ifndef ASTARIX_ARGPARSE_H
#define ASTARIX_ARGPARSE_H

#include <argp.h>
#include <cassert>
#include <cstring>
#include <string>

#include "utils.h"

// based on http://www.gnu.org/software/libc/manual/html_node/Argp-Example-3.html#Argp-Example-3

/* The options we understand. */
static struct argp_option options[] = {
    { "graph",          'g', "GRAPH",         0,  "Input graph (.gfa)" },
    { "query",          'q', "QUERY",         0,  "Input queries/reads (.fq, .fastq)" },
    { "outdir",         'o', "OUTDIR",        0,  "Output directory" },
    { "tree_depth",     'D', "TREE_DEPTH",    0,  "Suffix tree depth" },
    { "algorithm",      'a', "{dijkstra, astar-prefix}", 0, "Shortest path algorithm" },
    { "greedy_match",	'f', "GREEDY_MATCH",  0,  "Proceed greedily forward if there is a unique matching outgoing edge" },
    { "astar_len_cap",  'd', "A*_PREFIX_CAP", 0,  "The upcoming sequence length cap for the A* heuristic" },
    { "astar_cost_cap", 'c', "A*_COST_CAP",   0,  "The maximum prefix cost for the A* heuristic" },
    { "astar_equivalence_classes",
						'e', "A*_EQ_CLASSES", 0, "Whether to partition all nodes to equivalence classes in order not to reuse the heuristic" },
//    { "astar_lazy",		'L', "A*_LAZY",       0,  "Compute A* costs lazily during mapping" },
    { "match",			'M', "MATCH_COST",   0,  "Match penalty" },
    { "subst",	'S', "SUBST_COST",   0,  "Substitution penalty" },
    { "gap",			'G', "GAP_COST",     0,  "Gap (Insertion or Deletion) penalty" },
    { "threads",		't', "THREADS",      0,  "Number of threads (default=1)" },
    { "verbose",		'v', "THREADS",      0,  "Verbosity (default=silent=0, info=1, debug=2)" },
    { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments {
    char *command;
    char *graph_file;
    char *query_file;
    char *output_dir;
    char *algorithm;

	astarix::EditCosts costs;

	bool greedy_match;
    int tree_depth;

	int AStarLengthCap;
	double AStarCostCap;
	bool AStarNodeEqivClasses;

	int threads;
	int verbose;
};

/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    /* Get the input argument from argp_parse, which we
       know is a pointer to our arguments structure. */
    struct arguments *arguments = (struct arguments *)(state->input);
  
    switch (key) {
        case 'g':
            arguments->graph_file = arg;
            break;
        case 'q':
            arguments->query_file = arg;
            break;
        case 'D':
            arguments->tree_depth = std::stoi(arg);
            break;
        case 'a':
            assert(std::strcmp(arg, "dijkstra") == 0 || std::strcmp(arg, "astar-prefix") == 0);
            arguments->algorithm = arg;
            break;
		case 'f':
			arguments->greedy_match = (bool)std::stod(arg);
			break;
        case 'd':
			assert(std::stoi(arg) >= 0);
            arguments->AStarLengthCap = std::stoi(arg);
			break;
        case 'c':
			assert(std::stod(arg) >= 0.0);
			arguments->AStarCostCap = std::stod(arg);
			break;
		case 'e':
			arguments->AStarNodeEqivClasses = (bool)std::stod(arg);
			break;
        case 'o':
            arguments->output_dir = arg;
            break;
        case 'M':
            arguments->costs.match = std::stod(arg);
            break;
        case 'S':
            arguments->costs.subst = std::stod(arg);
            break;
        case 'G':
            arguments->costs.ins = std::stod(arg);
            arguments->costs.del = std::stod(arg);
            break;
        case 't':
            arguments->threads = std::stod(arg);
            break;
        case 'v':
            arguments->verbose = std::stod(arg);
            break;
        case ARGP_KEY_ARG:
            // Too many arguments.
            if (state->arg_num >= 3)
                argp_usage(state);
            assert(std::strcmp(arg, "align-optimal") == 0);
			arguments->command = arg;
            break;
  
        case ARGP_KEY_END:
            if (state->arg_num < 1)
                argp_usage(state);
            break;
  
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

arguments read_args(int argc, char **argv);

#endif
