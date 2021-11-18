#include <string>

#include "argparse.h"

struct argp_option options[] = {
    { "graph",          'g', "GRAPH",         0,  "Input graph (.gfa)" },
    { "query",          'q', "QUERY",         0,  "Input queries/reads (.fq, .fastq)" },
    { "outdir",         'o', "OUTDIR",        0,  "Output directory" },
    { "tree_depth",     'D', "TREE_DEPTH",    0,  "Suffix tree depth" },
    { "fixed_trie_depth",1001, "FIXED_TRIE_DEPTH",    0,  "Some leafs depth can be less than tree_depth (variable=0, fixed=1)" },
    { "algorithm",      'a', "{dijkstra, astar-prefix, astar-seeds}", 0, "Shortest path algorithm" },
    { "greedy_match",   'f', "GREEDY_MATCH",  0,  "Proceed greedily forward if there is a unique matching outgoing edge" },
    { "prefix_len_cap",  'd', "A*_PREFIX_CAP", 0,  "The upcoming sequence length cap for the A* heuristic" },
    { "prefix_cost_cap", 'c', "A*_COST_CAP",   0,  "The maximum prefix cost for the A* heuristic" },
    { "prefix_equivalence_classes",
                        'e', "A*_EQ_CLASSES", 0, "Whether to partition all nodes to equivalence classes in order not to reuse the heuristic" },
    { "seeds_len",  					2001, "A*_SEED_LEN", 0,  "The length of the A* seeds." },
    { "seeds_skip_near_crumbs",  		2008, "{0,1}", 0,  "" },
    { "match",          'M', "MATCH_COST",   0,  "Match penalty [0]" },
    { "subst",          'S', "SUBST_COST",   0,  "Substitution penalty [1]" },
    { "gap",            'G', "GAP_COST",     0,  "Gap (Insertion or Deletion) penalty [5]" },
    { "k_best_alignments", 'k', "TOP_K",     0,  "Output at most k optimal alignments per read [1]" },
    { "threads",        't', "THREADS",      0,  "Number of threads [1]" },
    { "verbose",        'v', "THREADS",      0,  "Verbosity (silent=0, info=1, debug=2), [0]" },
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
    args.graph_file            = "";
    args.query_file            = "";
    args.output_dir            = "";

    // Alignment parameters.
    args.costs                 = astarix::EditCosts(0, 1, 5, 5);
	args.k_best_alignments     = 1;               // top 5 optimal alignments

    // Performance parameters.
    args.algorithm             = (char *)"astar-prefix";
    args.tree_depth            = -1;              // auto mode
    args.fixed_trie_depth      = false;           // leafs can be shallower if `true`
    args.AStarLengthCap        = 5;
    args.AStarCostCap          = 5;
    args.threads               = 1;

    // Sound optimizations turned ON by default.
    args.greedy_match          = true;
    args.AStarNodeEqivClasses  = true;

    args.astar_seeds.seed_len              	= -1;
	args.astar_seeds.skip_near_crumbs		= true;

    args.verbose               = 0;
    args.command               = (char *)"align-optimal";
    
	
    /// ---- end of default values --------

    argp_parse(&argp, argc, argv, 0, 0, &args);

    if (args.graph_file == "") throw "Graph file not specified (-g).";
    if (args.query_file == "") throw "Query file not specified (-q).";

    if (!(args.costs.match   >= 0.0)) throw "EditCosts should be non-negative.";
    if (!(args.costs.subst   >= 0.0)) throw "EditCosts should be non-negative.";
    if (!(args.costs.ins     >= 0.0)) throw "EditCosts should be non-negative.";
    if (!(args.costs.del     >= 0.0)) throw "EditCosts should be non-negative.";
    if (!(args.costs.del     >= 0.0)) throw "EditCosts should be non-negative.";
    if (!(args.costs.match <= args.costs.subst)) throw "MatchCost should be not higher than SubstCost";
    if (!(args.costs.match <= args.costs.ins))   throw "MatchCost should be not higher than InsCost";
    if (!(args.costs.match <= args.costs.del))   throw "MatchCost should be not higher than DelCost";
    if (!(args.k_best_alignments >= 1)) throw "k_best_alignments should be at least 1.";

    if (!(args.threads >= 1)) throw "There should be a positive number of threads.";
    if (!(args.verbose >= 0)) throw "Verbosity should be non-negative.";

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
            if (std::strcmp(arguments->algorithm, "astar-prefix") != 0) throw "LengthCap only for astar-prefix.";
            if (!(std::stoi(arg) >= 0)) throw "AStarLengthCap should be non-negative.";
            arguments->AStarLengthCap = std::stoi(arg);
            break;
        case 'c':
            if (std::strcmp(arguments->algorithm, "astar-prefix") != 0) throw "CostCap only for astar-prefix.";
            if (!(std::stod(arg) >= 0.0)) throw "AStarCostCap should be non-negative.";
            arguments->AStarCostCap = std::stod(arg);
            break;
        case 'e':
            if (std::strcmp(arguments->algorithm, "astar-prefix") != 0) throw "NodeEquivClasses only for astar-prefix.";
            arguments->AStarNodeEqivClasses = (bool)std::stod(arg);
            break;
        case 2001:
            if (std::strcmp(arguments->algorithm, "astar-seeds") != 0) throw "SeedLen only for astar-seeds.";
            if (!(std::stoi(arg) >= 5)) throw "AStarSeedLen should be at least 5.";
            arguments->astar_seeds.seed_len = std::stod(arg);
            break;
        case 2008:
            arguments->astar_seeds.skip_near_crumbs = (bool)std::stod(arg);
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
        case 'k':
            arguments->k_best_alignments = std::stod(arg);
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
            if (std::strcmp(arg, "align-optimal") != 0) throw "align-optimal is a necessary command.";
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
