#include <string>

#include "astar-seeds-approx.h"
#include "argparse.h"

struct argp_option options[] = {
    { "graph",          'g', "GRAPH",         0,  "Input graph (.gfa)" },
    { "query",          'q', "QUERY",         0,  "Input queries/reads (.fq, .fastq)" },
    { "outdir",         'o', "OUTDIR",        0,  "Output directory" },
    { "tree_depth",     'D', "TREE_DEPTH",    0,  "Suffix tree depth" },
    { "fixed_trie_depth",1001, "FIXED_TRIE_DEPTH",    0,  "Some leafs depth can be less than tree_depth (variable=0, fixed=1)" },
    { "algorithm",      'a', "{dijkstra, astar-prefix}", 0, "Shortest path algorithm" },
    { "greedy_match",   'f', "GREEDY_MATCH",  0,  "Proceed greedily forward if there is a unique matching outgoing edge" },
    { "astar_len_cap",  'd', "A*_PREFIX_CAP", 0,  "The upcoming sequence length cap for the A* heuristic" },
    { "astar_cost_cap", 'c', "A*_COST_CAP",   0,  "The maximum prefix cost for the A* heuristic" },
    { "astar_equivalence_classes",
                        'e', "A*_EQ_CLASSES", 0, "Whether to partition all nodes to equivalence classes in order not to reuse the heuristic" },
//    { "astar_lazy",       'L', "A*_LAZY",       0,  "Compute A* costs lazily during mapping" },
    { "astar_seeds_len",  2001, "A*_SEED_LEN", 0,  "The length of the A* seeds." },
    { "astar_seeds_max_errors",  2002, "A*_SEEDS_MAX_ERRORS", 0,  "The maximum number of errors to a seed that a match can have." },
    { "astar_seeds_max_indels",  2003, "A*_SEEDS_MAX_INDELS", 0,  "The maximum number of indels. Any read with higher score with be reported as unaligned." },
    { "astar_seeds_backwards_algo",  2004, "{dfs_for_linear, bfs, complex, topsort}", 0,  "Backwards algo for each seed match." },
    { "match",          'M', "MATCH_COST",   0,  "Match penalty" },
    { "subst",          'S', "SUBST_COST",   0,  "Substitution penalty" },
    { "gap",            'G', "GAP_COST",     0,  "Gap (Insertion or Deletion) penalty" },
    { "threads",        't', "THREADS",      0,  "Number of threads (default=1)" },
    { "verbose",        'v', "THREADS",      0,  "Verbosity (default=silent=0, info=1, debug=2)" },
    { 0 }
};

// documentation: --help and --version
const char *argp_program_version = "AStarix 1.0";
const char *argp_program_bug_address = "<ivanov@pesho.info>";
static char doc[] = "Optimal sequence-to-graph aligner based on A* shortest path.";

/* A description of the arguments we accept. */
static char args_doc[] = "align-optimal -g GRAPH.gfa -q READS.fq -o OUT_DIR/";

static struct argp argp = { options, parse_opt, args_doc, doc };

arguments read_args(int argc, char **argv) {
    struct arguments args;

    // IO
    args.graph_file            = NULL;
    args.query_file            = NULL;
    args.output_dir            = NULL;

    // Problem statement.
    args.command               = (char *)"align-optimal";
    args.costs                 = astarix::EditCosts(0, 1, 5, 5);

    // Optimization parameters.
    args.algorithm             = (char *)"astar-prefix";
    args.tree_depth            = -1;              // auto mode
    args.fixed_trie_depth      = false;           // leafs can be shallower if `true`
    args.AStarLengthCap        = 5;
    args.AStarCostCap          = 5;

    // Sound optimizations turned ON by default.
    args.greedy_match          = true;
    args.AStarNodeEqivClasses  = true;

    args.astar_seeds.seed_len        = 15;
    args.astar_seeds.max_seed_errors = 0;
    args.astar_seeds.max_indels      = 10;
    args.astar_seeds.backwards_algo  = astarix::AStarSeedsWithErrors::Args::backwards_algo_t::BFS;

    args.threads               = 1;
    args.verbose               = 0;
    
    // Unsound optimizations turned OFF by default.
    // None.


    /// ---- end of default values --------

    argp_parse(&argp, argc, argv, 0, 0, &args);

    assert(args.graph_file && "Graph file not specified (-g).");
    assert(args.query_file && "Query file not specified (-q).");

    assert(args.costs.match   >= 0.0 && "EditCosts should be non-negative.");
    assert(args.costs.subst   >= 0.0 && "EditCosts should be non-negative.");
    assert(args.costs.ins     >= 0.0 && "EditCosts should be non-negative.");
    assert(args.costs.del     >= 0.0 && "EditCosts should be non-negative.");
    assert(args.costs.del     >= 0.0 && "EditCosts should be non-negative.");
    assert(args.costs.match <= args.costs.subst && "Match should be cheaper than other operations.");
    assert(args.costs.match <= args.costs.ins   && "Match should be cheaper than other operations.");
    assert(args.costs.match <= args.costs.del   && "Match should be cheaper than other operations.");

    assert(args.threads >= 1   && "There should be a positive number of threads.");
    assert(args.verbose >= 0   && "Verbosity should be non-negative.");

    return args;
}

error_t parse_opt (int key, char *arg, struct argp_state *state) {
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
        case 1001:
            arguments->fixed_trie_depth = (bool)std::stod(arg);
            break;
        case 'a':
            //assert(std::strcmp(arg, "dijkstra") == 0 || std::strcmp(arg, "astar-prefix") == 0);
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
        case 2001:
            assert(std::stoi(arg) >= 5);
            arguments->astar_seeds.seed_len = std::stod(arg);
            break;
        case 2002:
            assert(std::stoi(arg) >= 0);
            arguments->astar_seeds.max_seed_errors = std::stod(arg);
            break;
        case 2003:
            assert(std::stoi(arg) >= 0);
            arguments->astar_seeds.max_indels = std::stod(arg);
            break;
        case 2004:
            arguments->astar_seeds.backwards_algo = astarix::AStarSeedsWithErrors::Args::name2backwards_algo(arg);
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
