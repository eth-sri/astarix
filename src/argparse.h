#pragma once

#include <argp.h>
#include <cassert>
#include <cstring>
#include <string>

#include "utils.h"

// based on http://www.gnu.org/software/libc/manual/html_node/Argp-Example-3.html#Argp-Example-3

// long name, short key, arg name, flags, doc, group (0 is default)
/* The options we understand. */
extern struct argp_option options[];

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
    bool fixed_trie_depth;

    // A*-prefix params
    int AStarLengthCap;
    double AStarCostCap;
    bool AStarNodeEqivClasses;

    // A*-seeds params
    struct AStarSeedsArgs {
        int seed_len;
        int max_seed_errors;
        int max_indels;
    } astar_seeds;

    int threads;
    int verbose;
};

error_t parse_opt (int key, char *arg, struct argp_state *state);
arguments read_args(int argc, char **argv);
