//#define NDEBUG

#include <algorithm>
#include <dirent.h>
#include <errno.h>
#include <stdexcept>
#include <map>
#include <thread>
#include <mutex>

// For handing interruptions with CTRL+C
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "align.h"
#include "argparse.h"
#include "concurrentqueue.h"
#include "graph.h"
#include "io.h"
#include "trie.h"

// A* heuristics
#include "dijkstra.h"
#include "astar-prefix.h"
#include "astar-seeds-approx.h"
#include "astar-seeds-exact.h"

using namespace std;
using namespace astarix;

bool started_aligning;
bool interrupted;

// plog
void init_logger(const char *log_fn, int verbose) {
    if (verbose > 0) {
        assert(log_fn);
        auto level = verbose == 1 ? plog::info : plog::verbose;

        //static plog::RollingFileAppender<plog::TxtFormatter> InfoFileAppender(log_fn, 0);
        //plog::init<1>(level, &InfoFileAppender);
        //plog::init(level).addAppender(plog::get<1>());

        plog::init(level, log_fn); //1000000);
    }
}

unique_ptr<AStarHeuristic> AStarHeuristicFactory(const graph_t &G, const arguments &args) {
    unique_ptr<AStarHeuristic> astar;
    string algo = args.algorithm;

    if (algo == "astar-prefix") {
        astar = make_unique<AStarPrefix>(G, args.costs, args.AStarLengthCap, args.AStarCostCap, args.AStarNodeEqivClasses);
    } else if (algo == "astar-seeds-exact") {
        //if (!args.fixed_trie_depth)
        //    throw invalid_argument("astar-seeds-exact algorithm can only be used with fixed_trie_depth flag on.");
        //if (args.astar_seeds_max_errors != 0) 
        //    throw invalid_argument("astar-seeds-exact needs astar_seeds_max_errors flag set to 0.");
		throw "Exact not up to date.";
        astar = make_unique<AStarSeedsExact>(G, args.costs, args.astar_seeds.seed_len, -1);
    } else if (algo == "astar-seeds") {
        //if (!args.fixed_trie_depth)
        //    throw invalid_argument("astar-seeds algorithm can only be used with fixed_trie_depth flag on.");

        AStarSeedsWithErrors::Args astar_seeds = args.astar_seeds;

        astar = make_unique<AStarSeedsWithErrors>(G, args.costs, astar_seeds);
    } else if (algo == "dijkstra") { 
        astar = make_unique<DijkstraDummy>();
    } else {
        cout << "No algorithm " << args.algorithm << endl;
        throw invalid_argument("Unknown algorithm.");
    }

    return astar;
}

arguments args;

void wrap_readmap(const read_t& r, string algo, string performance_file, Aligner *aligner, bool calc_mapping_cost,
        edge_path_t *best_path, double *pushed_rate_sum, double *popped_rate_sum, double *repeat_rate_sum, double *pushed_rate_max, double *popped_rate_max, double *repeat_rate_max, FILE *fout,
       Stats *global_stats) {
    std::vector<state_t> final_states;
    
    aligner->astar_before_every_alignment(&r);      // prepare read
    final_states = aligner->readmap(r, algo, args.k_best_alignments);  // align
	aligner->astar_after_every_alignment();         // return preparation to previous state

	if ((int)final_states.size() >= args.k_best_alignments) {
		LOG_DEBUG << r.s << " aligned >= " << args.k_best_alignments << " times.";
	}

	for (auto &final_state: final_states) {
		best_path->clear();
		aligner->get_best_path_to_state(final_state, best_path);

		if (!performance_file.empty()) {
			string precomp_str = "align";
			int L = r.len;
			char strand = '?';
			
			int start = best_path->back().to;   // meaningful only for fasta where the nodeid is equal to the fasta position
			assert(start > 0);
			assert(start <= 2*((int)aligner->graph().nodes()+5));

			assert (!aligner->graph().node_in_trie(start));
			if (aligner->graph().node_in_reverse(start)) {
				start = aligner->graph().reverse2streight(start);
				//start -= L;
				strand = '-';
			} else {
				start -= L;
				strand = '+';
			}

	//		if (start > (int)aligner->graph().orig_nodes)
	//			start = -10000;

			double pushed_rate = (double)aligner->stats.pushed.get() / L;
			double popped_rate = (double)aligner->stats.popped.get() / L;
			double repeat_rate = (double)aligner->stats.repeated_visits.get() / aligner->stats.pushed.get();
			*pushed_rate_sum += pushed_rate;
			*popped_rate_sum += popped_rate;
			*repeat_rate_sum += repeat_rate;
			*pushed_rate_max = max(*pushed_rate_max, pushed_rate);
			*popped_rate_max = max(*popped_rate_max, popped_rate);
			*repeat_rate_max = max(*repeat_rate_max, repeat_rate);

			*global_stats += aligner->stats;
            int crumbs = aligner->astar->states_with_crumbs.get();

			char line[100000];
			line[0] = 0;
			sprintf(line,
					"%8s\t%3d\t%8s\t"
					"%8s\t%15s\t%8lf\t"
					"%3d\t%10s\t%10s\t"
					"%d\t%6d\t%c\t%6lf\t"
					"%6lf\t%4lf\t%8lf\t"
					"%8lf\t%d\t%d\t"
                    "%d\n",
					args.graph_file, (int)aligner->graph().nodes(), algo.c_str(),
					precomp_str.c_str(), r.comment.c_str(), 0.0,
					L, r.s.c_str(), spell(*best_path).c_str(),
					int(aligner->stats.align_status.cost.get()), start, strand, pushed_rate,
					popped_rate, repeat_rate, aligner->stats.t.total.get_sec(),
					aligner->stats.t.astar.get_sec(), aligner->stats.align_status.unique.get(), aligner->stats.explored_states.get(),
                    crumbs);
            fprintf(fout, "%s", line);
            fflush(fout);
		}
	}
}

void read_queries(const char *query_file, vector<read_t> *R) {
    read_t r;
    ifstream query_in(query_file);

    for (int i=0; read_query(query_in, query_file, &r); i++) {
        R->push_back(r);
        LOG_INFO_IF(i<5) << "read " << r.comment << ": " << r.s;
    }
    LOG_INFO_IF(R->size() >= 5) << " ...and more reads.";

    query_in.close();

    LOG_INFO << R->size() << " reads loaded.";
}

int size_sum(const vector<read_t> &R) {
    int sum = 0;
    for (const auto &r: R)
        sum += r.size();
    return sum;
}

double avg_error_rate(const vector<read_t> &R) {
    int letters = 0;
    double sump = 0.0;
    for (const auto &r: R) {
        for (size_t i=1; i<r.phreds.size(); i++) {
            int q = int(r.phreds[i])-33;
            double p = pow(10.0, -q/10.0);
            sump += p;
        }
        letters += r.size();
    }
    return sump/letters;
}

void auto_params(const graph_t &G, const vector<read_t> &R, arguments *args) {
    if (args->tree_depth == -1) {
        args->tree_depth = floor(log(G.nodes()) / log(4.0));
    }
    //if (args->astar_seeds.seed_len == -1) {
	//	args->astar_seeds.seed_len = args->tree_depth;
	//}
    if (args->tree_depth <= 0)
        throw "Trie depth should be >0.";
}

void print_tsv(map<string, string> dict, ostream &out) {
    for (auto const &it: dict)
        out << it.first << "\t";
    out << endl;
    for (auto const &it: dict) {
        assert(!it.second.empty());
        out << it.second << "\t";
    }
    out << endl;
}

//void print_hist(Aligner aligner, string hist_file) {
//    ofstream out(hist_file);
//    auto &all_counters = aligner.all_stats_;
//
//    int max_size = 0;
//    for (auto const &counters: all_counters) 
//        max_size = max(max_size, (int)counters.second.pushed_hist.size());
//
//    out << "read";
//    for (int i=0; i<max_size; i++)
//        out << "\t" << i;
//    out << "\n";
//
//    for (auto const &counters: all_counters) {
//        out << counters.first;
//        int i=0;
//        for (auto const &x: counters.second.pushed_hist) {
//            out << "\t" << x.get();
//            ++i;
//        }
//        for (; i<max_size; i++)
//            out << "\t" << 0;
//        out << "\n";
//    }
//    out << endl;
//    out.close();
//}

typedef map<string, string> dict_t;

struct Measurers {
    struct TimeAndMemory {
        Timer t;
        MemoryMeasurer m;

        void start() {
            t.start();
            m.start();
        }

        void stop() {
            t.stop();
            m.stop();
        }
    };

    TimeAndMemory total, construct_trie, read_graph, read_queries, align, precompute;

    void extract_to_dict(dict_t *dict) {
        (*dict)["total_sec"] = to_string(total.t.get_sec());
        (*dict)["total_gb"] = to_string(total.m.get_gb());
        (*dict)["construct_trie_sec"] = to_string(construct_trie.t.get_sec());
        (*dict)["construct_trie_gb"] = to_string(construct_trie.m.get_gb());
        (*dict)["read_graph_sec"] = to_string(read_graph.t.get_sec());
        (*dict)["read_graph_gb"] = to_string(read_graph.m.get_gb());
        (*dict)["read_queries_sec"] = to_string(read_queries.t.get_sec());
        (*dict)["read_queries_gb"] = to_string(read_queries.m.get_gb());
        (*dict)["align_sec"] = to_string(align.t.get_sec());
        (*dict)["align_gb"] = to_string(align.m.get_gb());
        (*dict)["precompute_sec"] = to_string(precompute.t.get_sec());
        (*dict)["precompute_gb"] = to_string(precompute.m.get_gb());
    }
};

void extract_args_to_dict(const arguments &args, dict_t *dict) {
    // io
    (*dict)["graph_file"] = args.graph_file;
    (*dict)["query_file"] = args.query_file;
    (*dict)["algorithm"] = args.algorithm;

    (*dict)["algorithm"] = args.algorithm;

    // edit costs
    (*dict)["cost_match"] = to_string(args.costs.match);
    (*dict)["cost_subst"] = to_string(args.costs.subst);
    (*dict)["cost_ins"] = to_string(args.costs.ins);
    (*dict)["cost_del"] = to_string(args.costs.del);

    // optimizations
    (*dict)["greedy_math"] = to_string(args.greedy_match);
    (*dict)["tree_depth"] = to_string(args.tree_depth);
    (*dict)["AStarLengthCap"] = to_string(args.AStarLengthCap);
    (*dict)["AStarCostCap"] = to_string(args.AStarCostCap);
    (*dict)["AStarNodeEqivClasses"] = to_string(args.AStarNodeEqivClasses);

    // perf
    (*dict)["threads"] = to_string(args.threads);
}

int exec(int argc, char **argv) {
    cout << "----" << endl;
#ifdef NDEBUG
    cout << "Assert() checks:       OFF" << endl;
#else
    cout << "Assert() checks:       ON" << endl;
#endif
    args = read_args(argc, argv);
    std::ios_base::sync_with_stdio(false);

    cout << "        verbose:        " << args.verbose << endl;
    cout << "----" << endl;

    string performance_file, info_log_file, stats_file, hist_file;

    string output_dir = args.output_dir;
    if (!output_dir.empty()) {
        assure_dir_exists(output_dir.c_str());
        performance_file = output_dir + "/alignments.tsv";
        info_log_file = output_dir + "/info.log";
        stats_file = output_dir + "/stats.log";
        hist_file = output_dir + "/hist.log";
    //if (!output_dir.empty())
        init_logger(info_log_file.c_str(), args.verbose);
    }

    Measurers T;
    dict_t stats;   // string key -> string value

    T.total.start();
    auto start_wt = std::chrono::high_resolution_clock::now();

    LOG_INFO << " ------------------------------ ";
    LOG_INFO << "Starting " << to_str(argc, argv);
    LOG_INFO << " ------------------------------ ";
    LOG_INFO << " sizeof(edge_t) = " << sizeof(edge_t);
    //LOG_DEBUG << "memory " << MemoryMeasurer::get_mem_gb();

    if (!performance_file.empty()) {
        FILE *fout = fopen(performance_file.c_str(), "w");
        fprintf(fout, "ref\trefsize\talgo\t"
                "operation\treadname\tmemory\t"
                "len\tread\tspell\t"
                "cost\tstart\tstrand\tpushed\t"
                "popped\trepeat_rate\tt(map)\t"
                "t(astar)\tunique_best\texplored_states\t"
                "crumbs\n");
        fclose(fout);
    }

    graph_t G;
    vector<read_t> R;
    //clock_t start;

    cout << "Loading reference graph... " << flush;
    T.read_graph.start();
    read_graph(&G, args.graph_file, output_dir);
	G.add_reverse_complement();
    cout << "Added reverse complement... " << flush;
    T.read_graph.stop();
    cout << "done in " << T.read_graph.t.get_sec() << "s."  << endl << flush;

    cout << "Loading queries... " << flush;
    T.read_queries.start();
    read_queries(args.query_file, &R);
    T.read_queries.stop();
    cout << "done in " << T.read_queries.t.get_sec() << "s." << endl << flush;

    auto_params(G, R, &args);

    cout << "Contructing trie... " << flush;
    T.construct_trie.start();
    add_tree(&G, args.tree_depth, args.fixed_trie_depth);
    T.construct_trie.stop();
    cout << "done in " << T.construct_trie.t.get_sec() << "s." << endl << flush;

    cout << "Initializing A* heuristic... " << flush;
    T.precompute.start();
    unique_ptr<AStarHeuristic> astar = AStarHeuristicFactory(G, args);  // TODO: resolve races
    T.precompute.stop();
    cout << "done in " << T.precompute.t.get_sec() << "s." << endl << flush;

    cost_t max_align_cost = 10000000; // args.astar_seeds.max_indels * std::min(args.costs.ins, args.costs.del);  // TODO: set with argument; not used
    AlignParams align_params(args.costs, args.greedy_match, max_align_cost);
    string algo = string(args.algorithm);

    assert(G.has_supersource());
    LOG_INFO << "Mapping init with graph with n=" << G.V.size() << " and m=" << G.E.size();
    align_params.print();

    std::ostream &out = cout;
    {
        out.setf(ios::fixed, ios::floatfield);
        out.precision(2);
        out << endl;
        out << " == General parameters and optimizations == "                                     << endl;
        out << "             Alignment algo: " << args.algorithm                                        << endl;
        out << "                 Edit costs: " << int(args.costs.match) << ", " << int(args.costs.subst) << ", "
                                    << int(args.costs.ins) << ", " << int(args.costs.del) << " (match, subst, ins, del)" << endl;
        out << "              Greedy match?: " << bool2str(args.greedy_match)                           << endl;
        out << "                    Threads: " << args.threads                                          << endl;
        out << endl;
        out << " == A* parameters =="                                                               << endl;
        astar->print_params(out);
        out << endl;
        out << " == Data =="                                                                      << endl;
        // Note: the trie is built on top of the **doubled** original graph (incl. reverse).
        out << "         Original reference: " << G.orig_nodes << " nodes, " << G.orig_edges << " edges"<< endl;
        out << "                       Trie: " << G.trie_nodes << " nodes, " << G.trie_edges << " edges, "
                                                << "depth" << (args.fixed_trie_depth ? "=" : "<=") << args.tree_depth                         << endl;
        out << "  Reference+ReverseRef+Trie: " << G.nodes() << " nodes, " << G.edges() << " edges, "
                                                << "density: " << (G.edges() / 2) / (G.nodes() / 2 * G.nodes() / 2) << endl;
        out << "                      Reads: " << R.size() << " x " << size_sum(R)/R.size() << "bp, "
                "coverage: " << 1.0 * size_sum(R) / ((G.edges() - G.trie_edges) / 2)<< "x" << endl; // the graph also includes reverse edges
        out << "            Avg phred value: " << 100.0*avg_error_rate(R) << "%" << endl;
        out << endl;

        stats["orig_graph_nodes"] = to_string(G.orig_nodes);
        stats["orig_graph_edges"] = to_string(G.orig_edges);
        stats["trie_nodes"] = to_string(G.trie_nodes); 
        stats["trie_edges"] = to_string(G.trie_edges);
        stats["total_nodes"] = to_string(G.nodes()); 
        stats["total_edges"] = to_string(G.edges());
    }

    double pushed_rate_sum(0.0), pushed_rate_max(0.0);
    double popped_rate_sum(0.0), popped_rate_max(0.0);
    double repeat_rate_sum(0.0), repeat_rate_max(0.0);

    T.align.start();
    auto start_align_wt = std::chrono::high_resolution_clock::now();

    std::mutex timer_m;
    atomic_int popped_trie_total(0), popped_ref_total(0);
    Stats global_stats;  // TODO: take care of races

    assert(G.E.size() == G.E_rev.size());  // TODO: remove
    assert(G.V.size() == G.V_rev.size());  // TODO: remove

    started_aligning = true;  // used for interruptions

    //cout << "Aligning..." << flush;
    bool calc_mapping_cost = false;
    if (args.threads == 1) {
        FILE *fout = fopen(performance_file.c_str(), "a");
        Aligner aligner(G, align_params, astar.get());
        for (size_t i=0; i<R.size(); i++) {
            if (interrupted)
                break;

            wrap_readmap(R[i], algo, performance_file, &aligner, calc_mapping_cost,
                    &R[i].edge_path, &pushed_rate_sum, &popped_rate_sum, &repeat_rate_sum, &pushed_rate_max, &popped_rate_max, &repeat_rate_max, fout, &global_stats);
            //global_stats.t += aligner.read_timers;

            popped_trie_total.fetch_add( aligner.stats.popped_trie.get() );  
            popped_ref_total.fetch_add( aligner.stats.popped_ref.get() );

            //if (i % (R.size() / 10) == 0) {
            //  cout << "A*-memoization at " << 100.0 * i / R.size() << "% of the reads aligned"
            //  << ", entries: " << astar.entries() << ", "
            //  << 100.0*b2gb(astar.table_mem_bytes_lower()) / MemoryMeasurer::get_mem_gb() << "%-"
            //  << 100.0*b2gb(astar.table_mem_bytes_upper()) / MemoryMeasurer::get_mem_gb() << "%" << endl;
            //}
        }
        fclose(fout);

        //print_hist(aligner, hist_file);
    } else {
        // TODO: handle interruptions using `interrupted` variable
        throw "Only 1 thread is currently supported.";  // TODO

        moodycamel::ConcurrentQueue<string> profileQueue { 50, (size_t)args.threads, (size_t)args.threads };
        std::vector<thread> threads(args.threads);
        std::atomic<bool> allThreadsDone { false };

        int bucket_sz = R.size() / args.threads;
        for (int t = 0; t < args.threads; ++t) {
            threads[t] = thread([&, t]() {
                int from = t*bucket_sz;
                int to = (t < args.threads-1) ? (t+1)*bucket_sz : R.size();
                unique_ptr<AStarHeuristic> astar_local = AStarHeuristicFactory(G, args);
                Aligner aligner(G, align_params, astar_local.get());
                LOG_INFO << "thread " << t << " for reads [" << from << ", " << to << ")";
                for (int i=from; i<to; i++) {
                    wrap_readmap(R[i], algo, performance_file, &aligner, calc_mapping_cost,
                            &R[i].edge_path, &pushed_rate_sum, &popped_rate_sum, &repeat_rate_sum, &pushed_rate_max, &popped_rate_max, &repeat_rate_max, NULL, &global_stats);
                    //profileQueue.enqueue(string(line));
                    //{
                    //    // TODO: merge the astar_local to astar stats
                    //    timer_m.lock();
            
                    //    timer_m.unlock();
                    //    if (t == 0) {
                    //        popped_trie_total.fetch_add( aligner.stats.popped_trie.get() );  
                    //        popped_ref_total.fetch_add( aligner.stats.popped_ref.get() );
                    //    }
                    //}
                }
            });
        }

        FILE *fout = fopen(performance_file.c_str(), "a");

        thread profileWriter([&]() {
            string line;
            while (true) {
                int elems = profileQueue.try_dequeue(line);
                if (!elems) {
                    if (allThreadsDone)
                        break;
                    else
                        std::this_thread::sleep_for(std::chrono::milliseconds(10));
                } else {
                    fprintf(fout, "%s", line.c_str());
                }
            }
        });

        for (int t = 0; t != args.threads; ++t) {
            threads[t].join();
        }

        allThreadsDone = true;

        profileWriter.join();
        fclose(fout);
    }
    T.align.stop();
    T.total.stop();
    auto end_align_wt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> align_wt = end_align_wt - start_align_wt;
    std::chrono::duration<double> total_wt = end_align_wt - start_wt;
    //cout << "done in " << T.align.t.get_sec() << "s." << endl << flush;

    {
        double total_mem = MemoryMeasurer::get_mem_gb();
        double align_cpu_time = T.align.t.get_sec();
        out << " == Aligning statistics =="                                                     << endl;
        out << "        Explored rate (avg): " << 1.0*global_stats.explored_states.get() / R.size() / R[0].len << " states/read_bp" << endl;
        out << "     Pushed rate (avg, max): " << pushed_rate_sum / R.size() << ", " << pushed_rate_max << "    [states/bp] (states normalized by query length)" << endl;
        out << "     Popped rate (avg, max): " << popped_rate_sum / R.size() << ", " << popped_rate_max << endl;
        out << "             Average popped: " << 1.0 * popped_trie_total.load() / (R.size()/args.threads)
                                            << " from trie (" << 100.0*popped_trie_total.load()/(popped_trie_total.load() + popped_ref_total.load()) << "%) vs "
                                            << 1.0 * popped_ref_total.load() / (R.size()/args.threads) << " from ref"  << " (per read)" << endl;
        out << "Total cost of aligned reads: " << global_stats.align_status.cost.get() << ", " << 1.*global_stats.align_status.cost.get()/global_stats.align_status.aligned() << " per read, " 
            << 100.0*global_stats.align_status.cost.get()/size_sum(R) << "% per letter" << endl;
#ifndef NDEBUG
        out << "      Repeated states (avg): " << 1.0*global_stats.repeated_visits.get() / R.size() / R[0].len << " states/read_bp" << endl;
#endif
        out << "                 Alignments: " 
                                                << global_stats.align_status.unique.get()    << " unique ("    << 100.*global_stats.align_status.unique.get()/R.size()    << "%), "
                                                << global_stats.align_status.ambiguous.get() << " ambiguous (" << 100.*global_stats.align_status.ambiguous.get()/R.size() << "%) and "
                                                << global_stats.align_status.overcost.get()  << " overcost ("  << 100.*global_stats.align_status.overcost.get()/R.size()  << "%)" << endl;
        out << endl;
        out << " == Heuristic stats (" << args.algorithm << ") ==" << std::endl;
        astar->print_stats(out);
        out << endl;
        out << " == Performance =="                                                             << endl;
        out << "    Memory: " << "                   measured | estimated"                                  << endl;
        out << "                   total: " << total_mem << "gb, 100% | -"      << endl;
        out << "               reference: " << T.read_graph.m.get_gb() << "gb, " << 100.0*T.read_graph.m.get_gb() / total_mem << "% | " << 100.0*b2gb(G.reference_mem_bytes()) / total_mem << "%" << endl;
        out << "                   reads: " << T.read_queries.m.get_gb() << "gb, " << 100.0*T.read_queries.m.get_gb() / total_mem << "% | " << 100.0*b2gb(R.size() * R.front().size()) / total_mem << "%" << endl;
        out << "                    trie: " << T.construct_trie.m.get_gb() << "gb, " << 100.0*T.construct_trie.m.get_gb() / total_mem << "% | " << 100.0*b2gb(G.trie_mem_bytes()) / total_mem << "%" << endl;
        out << "     equiv. classes opt.: " << T.precompute.m.get_gb() << "gb, " << 100.0*T.precompute.m.get_gb() / total_mem << "%" << endl;
        out << "          A*-memoization: " << T.align.m.get_gb() << "gb, " << 100.0*T.align.m.get_gb() / total_mem << endl;
        out << endl;
        out << "   Total wall runtime:    " << total_wt.count() << "s"                  << endl;
        out << "       reference loading: " << T.read_graph.t.get_sec() << "s"          << endl;
        out << "         queries loading: " << T.read_queries.t.get_sec() << "s"        << endl;
        out << "          construct trie: " << T.construct_trie.t.get_sec() << "s"      << endl;
        out << "              precompute: " << T.precompute.t.get_sec() << "s"          << endl;
        out << "       align (wall time): " << align_wt.count() << "s = "
                                            << R.size() / align_wt.count() << " reads/s = "
                                            << size_sum(R) / 1000.0 / align_wt.count() << " Kbp/s"    << endl; 
        out << endl;
        out << "    Total align cpu time: " << align_cpu_time << "s = "
                                            << R.size() / align_cpu_time << " reads/s = "
                                            << size_sum(R) / 1000.0 / align_cpu_time << " Kbp/s"    << endl; 
        out << "     |          Preprocessing: " << 100.0 * global_stats.t.astar_prepare_reads.get_sec() / align_cpu_time << "%" << endl;
        out << "     |               A* query: " << 100.0 * global_stats.t.astar.get_sec() / align_cpu_time << "%"   << endl;
        out << "     |                  queue: " << 100.0 * global_stats.t.queue.get_sec() / align_cpu_time << "%" << endl;
        out << "     |                  dicts: " << 100.0 * global_stats.t.dicts.get_sec() / align_cpu_time << "%" << endl;
        out << "     |           greedy_match: " << 100.0 * global_stats.t.ff.get_sec() / align_cpu_time << "%" << endl;
        out << "     |                  other: " << 100.0 - 100.0 * (global_stats.t.astar_prepare_reads.get_sec()
                                                                + global_stats.t.astar.get_sec()
                                                                + global_stats.t.queue.get_sec()
                                                                + global_stats.t.dicts.get_sec()
                                                                + global_stats.t.ff.get_sec()) / align_cpu_time << "%" << endl;
        out << " DONE" << endl;
        out << endl;

        assert( global_stats.align_status.total() >= (int)R.size() );
    }

    extract_args_to_dict(args, &stats);
    T.extract_to_dict(&stats);

    {
        ofstream tsv(stats_file);
        print_tsv(stats, tsv);
        tsv.close();
    }


    return 0;
}

void my_handler(int s) {
    printf("Caught signal %d\n", s);
    interrupted = true;

    if (!started_aligning)
        exit(1); 
}

int main(int argc, char **argv) {
    {
        // Catching CTRL+C
        started_aligning = false;
        interrupted = false;
        struct sigaction sigIntHandler;

        sigIntHandler.sa_handler = my_handler;
        sigemptyset(&sigIntHandler.sa_mask);
        sigIntHandler.sa_flags = 0;

        sigaction(SIGINT, &sigIntHandler, NULL);
    }

    try {
        return exec(argc, argv);
    } catch (const std::string& ex) {
        std::cout << std::endl;
        std::cout << "Caught exception: " << ex << std::endl;
        LOG_FATAL << ex;
		return 1;
    } catch (const char *ex) {
        std::cout << std::endl;
        std::cout << "Caught exception: " << ex << std::endl;
        LOG_FATAL << ex;
		return 1;
    }
}

