<img width="100" alt="AStarix" align="left" src="https://www.sri.inf.ethz.ch/assets/systems/astarix.png"><br/>
<a href="https://www.sri.inf.ethz.ch/"><img width="100" alt="SRILAB" align="right" src="http://safeai.ethz.ch/img/sri-logo.svg"></a><br/>

---

![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)

**AStarix** is a sequence-to-graph aligner based on A* shortest path algorithm.
It supports general graphs and finds alignments that are optimal according to edit-distance with non-negative weights.

## Installation

Alternatively to using Docker, you can install `AStarix` directly on your system.

### Prerequisites

In order to compile `AStarix`, you will need to first install:

* [argp](https://www.gnu.org/software/libc/manual/html_node/Argp.html) &ndash;
  argument parsing library
* (optional) [pandas](https://pandas.pydata.org/) &ndash; dataframe library used for testing and evalutaions

T hird-party libraries are located in the `/ext` directory and their own licenses apply.

### Compile and test

```bash
make
make test
```

# Example run

This is an example run of `AStarix` on a toy data: aligning 100 simulated Illumina reads to a linear graph from the first 10000bp of Escherichia coli. You can expect a similar summary printed to standard output:

```
$ make test
release/astarix align-optimal -a astar-seeds -t 1 -v 0 -g data/ecoli_head10000_linear/graph.gfa -q data/ecoli_head10000_linear/illumina.fq -o tmp/ecoli_head10000_linear/astar-seeds --fixed_trie_depth 1 --seeds_len 25
----
Assert() checks:       OFF
        verbose:        0
----
Loading reference graph... Added reverse complement... done in 0.001078s.
Loading queries... done in 0.000216s.
Contructing trie... done in 0.007017s.
Initializing A* heuristic... done in 1.6e-05s.

 == General parameters and optimizations ==
             Alignment algo: astar-seeds
                 Edit costs: 0, 1, 5, 5 (match, subst, ins, del)
              Greedy match?: true
                    Threads: 1

 == A* parameters ==
          seed length: 25 bp
     skip near crumbs: 1

 == Data ==
         Original reference: 9826 nodes, 9824 edges
                       Trie: 5186 nodes, 24822 edges, depth=7
  Reference+ReverseRef+Trie: 24838 nodes, 44470 edges, density: 0

                  Trie: 5186 nodes, 24822 edges, depth=7
  Reference+ReverseRef+Trie: 24838 nodes, 44470 edges, density: 0
                      Reads: 100 x 99bp, coverage: 1.00774x
            Avg phred value: 0.13640%

 == Aligning statistics ==
        Explored rate (avg): 1.70707 states/read_bp
         States with crumbs: 0.02992%
            Explored states: 0.01737%
             Skipped states: 99.95271%
     Pushed rate (avg, max): 0.77130, 0.01400    [states/bp] (states normalized by query length)
     Popped rate (avg, max): 0.09310, 0.00140
             Average popped: 7.16000 from trie (76.90655%) vs 2.15000 from ref (per read)
Total cost of aligned reads: 15, 0.15000 per read, 0.15152% per letter
                 Alignments: 100 unique (100.00000%), 0 ambiguous (0.00000%) and 0 overcost (0.00000%)

 == Heuristic stats (astar-seeds) ==
        For all reads:
                            Seeds: 400 (4.00000 per read)
                     Seed matches: 385 (3.85000 per read, 0.96250 per seed)
               States with crumbs: 29105 [+0.00000% repeated], (291.05000 per read)
                  Heuristic (avg): 0.10000 of potential 4.00000

 == Performance ==
    Memory:                    measured | estimated
                   total: 0.00510gb, 100% | -
               reference: 0.00267gb, 52.39521% | 5.74505%
                   reads: 0.00000gb, 0.00000% | 0.18091%
                    trie: 0.00100gb, 19.53593% | 5.82224%
     equiv. classes opt.: 0.00000gb, 0.00000%
          A*-memoization: 0.00000gb, 0.00000

   Total wall runtime:    0.04086s
       reference loading: 0.00108s
         queries loading: 0.00022s
          construct trie: 0.00702s
              precompute: 0.00002s
       align (wall time): 0.03040s = 3289.04470 reads/s = 325.61543 Kbp/s

    Total align cpu time: 0.02998s = 3335.11206 reads/s = 330.17609 Kbp/s
     |          Preprocessing: 31.61686%
     |               A* query: 18.47652%
     |           greedy_match: 11.42609%
     |                  other: 38.48052%
 DONE
```

The directory `tmp/ecoli_head10000_linear/astar-seeds/` will contain a file with execution logs and a file with alignment statistics.
Short aggregated statistics are print to standard output (to redirect, you can add `>summary.txt`).


# Usage

`AStarix` currently works only in exact/optimal mode (specified using `align-optimal`). Currently supported formats are `.gfa` without overlapping nodes (for a graph reference) and `.fa`/`.fasta` (for a linear reference). The queries should be in `.fq`/`.fastq` format (the phred values are ignored).
Tested on Ubuntu 20.04.

```
$ astarix --help
----  
Usage: astarix [OPTION...] align-optimal -g GRAPH.gfa -q READS.fq -o OUT_DIR/
Optimal sequence-to-graph aligner based on A* shortest path.

  -a, --algorithm={dijkstra, astar-prefix, astar-seeds}
                             Shortest path algorithm
  -c, --prefix_cost_cap=A*_COST_CAP
                             The maximum prefix cost for the A* heuristic
  -d, --prefix_len_cap=A*_PREFIX_CAP
                             The upcoming sequence length cap for the A*
                             heuristic
  -D, --tree_depth=TREE_DEPTH   Suffix tree depth
  -e, --prefix_equivalence_classes=A*_EQ_CLASSES
                             Whether to partition all nodes to equivalence
                             classes in order not to reuse the heuristic
      --fixed_trie_depth=FIXED_TRIE_DEPTH
                             Some leafs depth can be less than tree_depth
                             (variable=0, fixed=1)
  -f, --greedy_match=GREEDY_MATCH
                             Proceed greedily forward if there is a unique
                             matching outgoing edge
  -g, --graph=GRAPH          Input graph (.gfa)
  -G, --gap=GAP_COST         Gap (Insertion or Deletion) penalty [5]
  -k, --k_best_alignments=TOP_K   Output at most k optimal alignments per read
                             [1]
  -M, --match=MATCH_COST     Match penalty [0]
  -o, --outdir=OUTDIR        Output directory
  -q, --query=QUERY          Input queries/reads (.fq, .fastq)
      --seeds_len=A*_SEED_LEN   The length of the A* seeds.
      --seeds_skip_near_crumbs={0,1}
  -S, --subst=SUBST_COST     Substitution penalty [1]
  -t, --threads=THREADS      Number of threads [1]
  -v, --verbose=THREADS      Verbosity (silent=0, info=1, debug=2), [0]
  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version
```

# Technical details

`AStarix` computes an A* heuristic function __h__ that directs the path search faster by anticipating the upcoming nucleotides to be aligned.
Once computed, a value of __h__ is memoized in a hash table ([parallel-hashmap](https://github.com/greg7mdp/parallel-hashmap)).

`AStarix` parallelizes the alignment of a set of reads.

## Citing

[1] [Trie and Prefix heuristic](https://www.biorxiv.org/content/10.1101/2020.01.22.915496v1) (2020) -- Introducing AStarix tool, the Trie index which allows scaling to bigger references, and Prefix heuristic which outperforms Dijkstra. The paper was presented on [RECOMB 2020](https://www.recomb2020.org/). Evaluations [here](https://github.com/eth-sri/astarix/tree/RECOMB2020_experiments/evals).

[2] [Seed heuristic preprint](https://www.biorxiv.org/content/10.1101/2021.11.05.467453v1) (2021) -- Introducing _Seed heuristic_ which speeds up aligning of short reads and allow scaling to long HiFi reads.

# Historical note
The development of the `AStarix` algorithm is only possible because of a
multitude of conceptual inventions, including:

* [Dijkstra's
  algorithm](https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm) (1959) &ndash;
  provably optimal (given non-negative edges), provably fast (polynomial worst
  case) greedy shortest path algorithm that considers first the nodes closest to
  the start
* [Levenshtein distance](https://en.wikipedia.org/wiki/Levenshtein_distance)
  (1965) &ndash; the variant of _edit distance_ that considers insertions, deletions
  and subtitutions
* [A* algorithm](https://en.wikipedia.org/wiki/A*_search_algorithm) (1968) &ndash; extension of Dijkstra's algorithm that also considers a heuristic
  estimate of the remaining path length; provably optimal (using an
  _admissibile_/_optimistic_ heuristic), not slower than Dijkstra in the worst
  case (assuming a _consistent_ heuristic), crucial for the result optimality,
  may enormously speed up Dijkstra's algorithm (in average case) when used with
  an informative heuristic
* [Suffix tree](https://en.wikipedia.org/wiki/Suffix_tree) (1973) &ndash; data
  structure used to directly search from all possible starting places of a string

Among the countless other indirect influences, most notable conceptual works include:
* [Alan Turing's early program proof
  work](https://fi.ort.edu.uy/innovaportal/file/20124/1/09-turing_checking_a_large_routine_earlyproof.pdf)
  (1949) 
* [Automata theory](https://en.wikipedia.org/wiki/Automata_theory)
* [Dynamic programming](https://en.wikipedia.org/wiki/Dynamic_programming)