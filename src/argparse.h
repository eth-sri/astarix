#pragma once

#include <argp.h>
#include <cassert>
#include <cstring>
#include <map>
#include <string>

#include "astar-seeds.h"
#include "utils.h"

// based on http://www.gnu.org/software/libc/manual/html_node/Argp-Example-3.html#Argp-Example-3

// long name, short key, arg name, flags, doc, group (0 is default)
/* The options we understand. */
extern struct argp_option options[];


/* Used by main to communicate with parse_opt. */
struct arguments {
    // Input
	std::string graph_file;
	std::string query_file;
	std::string output_dir;

	// Alignment parameters
    astarix::EditCosts costs;
	int k_best_alignments;

    // Performance params
    char *algorithm;
    bool greedy_match;
    int tree_depth;
    bool fixed_trie_depth;
    int threads;

    // A*-prefix params
    int AStarLengthCap;
    double AStarCostCap;
    bool AStarNodeEqivClasses;
    astarix::AStarSeedsWithErrors::Args astar_seeds;

    // Debug
    int verbose;
    char *command;
};

error_t parse_opt (int key, char *arg, struct argp_state *state);
arguments read_args(int argc, char **argv);
