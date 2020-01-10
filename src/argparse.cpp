#include "argparse.h"

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
	args.command			   = "align-optimal";
    args.costs                 = astarix::EditCosts(0, 1, 5, 5);

	// Optimization parameters.
	args.algorithm             = "astar-prefix";
	args.tree_depth            = -1;              // auto mode
	args.AStarLengthCap        = 5;
	args.AStarCostCap          = 5;

	// Sound optimizations turned ON by default.
	args.greedy_match          = true;
	args.AStarNodeEqivClasses  = true;

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
