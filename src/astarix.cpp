//#define NDEBUG

#include <algorithm>
#include <dirent.h>
#include <errno.h>
#include <stdexcept>
#include <map>
#include <thread>

#include "align.h"
#include "argparse.h"
#include "concurrentqueue.h"
#include "graph.h"
#include "io.h"
#include "trie.h"

// A* heuristics
#include "dijkstra.h"
#include "astar-prefix.h"
#include "astar-landmarks-approx.h"
#include "astar-landmarks-exact.h"

using namespace std;
using namespace astarix;

// plog
void init_logger(const char *log_fn, int verbose) {
    if (verbose > 0) {
        assert(log_fn);
        auto level = verbose == 1 ? plog::info : plog::verbose;

        //static plog::RollingFileAppender<plog::TxtFormatter> InfoFileAppender(log_fn, 0);
        //plog::init<1>(level, &InfoFileAppender);
        //plog::init(level).addAppender(plog::get<1>());

        plog::init(level, log_fn, 1000000);
    }
}

AStarHeuristic* AStarHeuristicFactory(const graph_t &G, const arguments &args) {
    AStarHeuristic* astar;
    string algo = args.algorithm;

    int shifts_allowed = 25;

    // TODO: add dijkstra
    if (algo == "astar-prefix") {
        astar = new AStarPrefix(G, args.costs, args.AStarLengthCap, args.AStarCostCap, args.AStarNodeEqivClasses);
    } else if (algo == "astar-landmarks-exact") {
        if (!args.fixed_trie_depth)
            throw invalid_argument("astar-landmarks-exact algorithm can only be used with fixed_trie_depth flag on.");
        if (args.astar_max_waymark_errors != 0) 
            throw invalid_argument("astar-landmarks-exact needs astar_max_waymark_errors flag set to 0.");
        astar = new AStarLandmarksExact(G, args.costs, args.astar_landmark_len, shifts_allowed);
    } else if (algo == "astar-landmarks") {
        if (!args.fixed_trie_depth)
            throw invalid_argument("astar-landmarks algorithm can only be used with fixed_trie_depth flag on.");
        astar = new AStarWaymarksWithErrors(G, args.costs, args.astar_landmark_len, args.astar_max_waymark_errors, shifts_allowed);
    } else if (algo == "dijkstra") { 
        astar = new DijkstraDummy();
    } else {
        cout << "No algorithm " << args.algorithm << endl;
        throw invalid_argument("Unknown algorithm.");
    }

    return astar;
}

arguments args;

state_t wrap_readmap(const read_t& r, string algo, string performance_file, Aligner *aligner, bool calc_mapping_cost,
        edge_path_t *path, double *pushed_rate_sum, double *popped_rate_sum, double *repeat_rate_sum, double *pushed_rate_max, double *popped_rate_max, double *repeat_rate_max, char *line,
       Counters *all_reads_counters) {
    state_t ans;
    
    // prepare read
    aligner->astar_before_every_alignment(&r);

    ans = aligner->readmap(r, algo, path);

    // return preparation to previous state
    aligner->astar_after_every_alignment();

    if (!performance_file.empty()) {
        //const auto &astar = aligner->get_astar();
        const auto &timers = aligner->read_timers;

        string precomp_str = "align";
        int L = r.len;
        int starts = -1;
        double pushed_rate = (double)aligner->read_counters.pushed.get() / L;
        double popped_rate = (double)aligner->read_counters.popped.get() / L;
        double repeat_rate = (double)aligner->_repeated_visits / aligner->read_counters.pushed.get();
        *pushed_rate_sum += pushed_rate;
        *popped_rate_sum += popped_rate;
        *repeat_rate_sum += repeat_rate;
        *pushed_rate_max = max(*pushed_rate_max, pushed_rate);
        *popped_rate_max = max(*popped_rate_max, popped_rate);
        *repeat_rate_max = max(*repeat_rate_max, repeat_rate);

        *all_reads_counters += aligner->read_counters;

        line[0] = 0;
        sprintf(line,
                "%8s\t%3d\t%8s\t"
                "%8s\t%15s\t%8lf\t"
                "%3d\t%10s\t%10s\t"
                "%d\t%6d\t%6lf\t"
                "%6lf\t%4lf\t%8lf\t"
                "%8lf\t%d\n",
                args.graph_file, (int)aligner->graph().nodes(), algo.c_str(),
                precomp_str.c_str(), r.comment.c_str(), 0.0,
                L, r.s.c_str(), spell(*path).c_str(),
                int(ans.cost), starts, pushed_rate,
                popped_rate, repeat_rate, timers.total.get_sec(),
                timers.astar.get_sec(), aligner->unique_best);
    }

    return ans;
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

void auto_params(const graph_t &G, const vector<read_t> &R, arguments *args) {
    if (args->tree_depth == -1) {
        args->tree_depth = floor(log(G.nodes()) / log(4.0));
    }
    assert(args->tree_depth > 0);
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

int main(int argc, char **argv) {
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

    string performance_file, info_log_file, stats_file;

    string output_dir = args.output_dir;
    if (!output_dir.empty()) {
        assure_dir_exists(args.output_dir);
        performance_file = output_dir + "/alignments.tsv";
        info_log_file = output_dir + "/info.log";
        stats_file = output_dir + "/stats.log";
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
                "cost\tstarts\tpushed\t"
                "popped\trepeat_rate\tt(map)\t"
                "t(astar)\tunique_best\n");
        fclose(fout);
    }

    graph_t G;
    vector<read_t> R;
    clock_t start;

    cout << "Loading reference graph... " << flush;
    T.read_graph.start();
    read_graph(&G, args.graph_file, output_dir);
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
    AStarHeuristic *astar = AStarHeuristicFactory(G, args);
    T.precompute.stop();
    cout << "done in " << T.precompute.t.get_sec() << "s." << endl << flush;

    AlignParams align_params(args.costs, args.greedy_match);
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
                                                << "depth: " << args.tree_depth                         << endl;
        out << "  Reference+ReverseRef+Trie: " << G.nodes() << " nodes, " << G.edges() << " edges, "
                                                << "density: " << (G.edges() / 2) / (G.nodes() / 2 * G.nodes() / 2) << endl;
        out << "                      Reads: " << R.size() << " x " << R.front().len << "bp, "
                "coverage: " << 1.0 * R.size() * (R.front().s.size()-1) / ((G.edges() - G.trie_edges) / 2)<< "x" << endl; // the graph also includes reverse edges
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

    AlignerTimers total_timers;
    int total_cost(0);
    std::mutex timer_m;
    Counter popped_trie_total, popped_ref_total;
    Counters all_reads_counters;
    std::atomic<bool> nonunique_best_alignments(0);

    assert(G.E.size() == G.E_rev.size());  // TODO: remove
    assert(G.V.size() == G.V_rev.size());  // TODO: remove

    cout << "Aligning..." << flush;
    bool calc_mapping_cost = false;
    if (args.threads == 1) {
        FILE *fout = fopen(performance_file.c_str(), "a");
        Aligner aligner(G, align_params, astar);
        for (size_t i=0; i<R.size(); i++) {
            char line[10000];
            state_t ans = wrap_readmap(R[i], algo, performance_file, &aligner, calc_mapping_cost,
                    &R[i].edge_path, &pushed_rate_sum, &popped_rate_sum, &repeat_rate_sum, &pushed_rate_max, &popped_rate_max, &repeat_rate_max, line, &all_reads_counters);
            fprintf(fout, "%s", line);
            total_timers += aligner.read_timers;
            total_cost += ans.cost;
            if (!aligner.unique_best)
                nonunique_best_alignments = nonunique_best_alignments+1;

            //if (i % (R.size() / 10) == 0) {
            //  cout << "A*-memoization at " << 100.0 * i / R.size() << "% of the reads aligned"
            //  << ", entries: " << astar.entries() << ", "
            //  << 100.0*b2gb(astar.table_mem_bytes_lower()) / MemoryMeasurer::get_mem_gb() << "%-"
            //  << 100.0*b2gb(astar.table_mem_bytes_upper()) / MemoryMeasurer::get_mem_gb() << "%" << endl;
            //}
        }
        fclose(fout);
    } else {
        moodycamel::ConcurrentQueue<string> profileQueue { 50, args.threads, args.threads };
        std::vector<thread> threads(args.threads);
        std::atomic<bool> allThreadsDone { false };

        int bucket_sz = R.size() / args.threads;
        for (int t = 0; t < args.threads; ++t) {
            threads[t] = thread([&, t]() {
                int from = t*bucket_sz;
                int to = (t < args.threads-1) ? (t+1)*bucket_sz : R.size();
                AStarHeuristic *astar_local = AStarHeuristicFactory(G, args);
                Aligner aligner(G, align_params, astar_local);
                LOG_INFO << "thread " << t << " for reads [" << from << ", " << to << ")";
                for (size_t i=from; i<to; i++) {
                    char line[10000];
                    state_t ans = wrap_readmap(R[i], algo, performance_file, &aligner, calc_mapping_cost,
                            &R[i].edge_path, &pushed_rate_sum, &popped_rate_sum, &repeat_rate_sum, &pushed_rate_max, &popped_rate_max, &repeat_rate_max, line, &all_reads_counters);
                    profileQueue.enqueue(string(line));
                    {
                        // TODO: merge the astar_local to astar stats
                        timer_m.lock();
                        total_cost += ans.cost;
                        total_timers += aligner.read_timers;
            
                        timer_m.unlock();
                        if (t == 0) {
                            popped_trie_total.inc( aligner.read_counters.popped_trie.get() );  
                            popped_ref_total.inc( aligner.read_counters.popped_ref.get() );
                        }
                        if (!aligner.unique_best)
                            nonunique_best_alignments = nonunique_best_alignments + 1;
                    }
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
    cout << "done in " << T.align.t.get_sec() << "s." << endl << flush;

    {
        double total_map_time = total_timers.total.get_sec();
        double total_mem = MemoryMeasurer::get_mem_gb();

        out << " == Aligning statistics =="                                                     << endl;
        out << "   Explored rate (avg, max): " << 1.0*all_reads_counters.explored_states.get() / R.size() / R[0].len << " states/read_bp" << endl;
        out << "    Pushed  rate (avg, max): " << pushed_rate_sum / R.size() << ", " << pushed_rate_max << "    [states/bp] (states normalized by query length)" << endl;
        out << "     Popped rate (avg, max): " << popped_rate_sum / R.size() << ", " << popped_rate_max << endl;
        out << "             Average popped: " << 1.0 * popped_trie_total.get() / (R.size()/args.threads) << " from trie vs "
                                            << 1.0 * popped_ref_total.get() / (R.size()/args.threads) << " from ref"
                                            << " (influenced by greedy match flag)" << endl;
        out << "                 Total cost: " << total_cost << endl;
        out << " Non-unique best alignments: " << nonunique_best_alignments << " out of " << R.size()   << endl;
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
        out << "                   align: " << align_wt.count() << "s (wall time) => "
                                            << R.size() / total_map_time << " reads/s <=> "
                                            << size_sum(R) / 1024.0 / total_map_time << " Kbp/s"    << endl; 
        out << endl;
        out << " Total align cpu time:    " << T.align.t.get_sec() << "s (cpu time)"    << endl;
        out << "                      A*: " << 100.0 * total_timers.astar.get_sec() / total_map_time   << "%, "   << endl;
        out << "        A* prepare_reads: " << 100.0 * total_timers.astar_prepare_reads.get_sec() / total_map_time    << "%, " << endl;
        out << "                     que: " << 100.0 * total_timers.queue.get_sec() / total_map_time    << "%, " << endl;
        out << "                   dicts: " << 100.0 * total_timers.dicts.get_sec() / total_map_time  << "%, " << endl;
        out << "            greedy_match: " << 100.0 * total_timers.ff.get_sec() / total_map_time  << "%" << endl;
        out << " DONE" << endl;
    }

    extract_args_to_dict(args, &stats);
    T.extract_to_dict(&stats);
    ofstream tsv(stats_file);
    print_tsv(stats, tsv);
    tsv.close();

    return 0;
}

