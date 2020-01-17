//#define NDEBUG

#include <algorithm>
#include <dirent.h>
#include <errno.h>
#include <map>
#include <thread>

#include "align.h"
#include "argparse.h"
#include "concurrentqueue.h"
#include "graph.h"
#include "io.h"
#include "trie.h"

using namespace std;
using namespace astarix;

arguments args;

// plog
void init_logger(const char *log_fn, int verbose) {
#ifndef NDEBUG
	if (log_fn && verbose > 0) {
		auto level = verbose == 1 ? plog::info : plog::debug;
		static plog::RollingFileAppender<plog::TxtFormatter> InfoFileAppender(log_fn, 0);
		plog::init<1>(level, &InfoFileAppender);
		
		plog::init(level).addAppender(plog::get<1>());
	}
#endif
}

state_t wrap_readmap(read_t& r, string algo, string performance_file, Aligner *aligner, bool calc_mapping_cost,
		edge_path_t *path, double *pushed_rate_sum, double *popped_rate_sum, double *repeat_rate_sum, double *pushed_rate_max, double *popped_rate_max, double *repeat_rate_max, char *line) {
	state_t ans;
	
	ans = aligner->readmap(r, algo, path);

	const auto &astar = aligner->get_astar();
	const auto &timers = aligner->read_timers;

	if (!performance_file.empty()) {
		string precomp_str = "align";
		string pref = r.s.substr(1, astar.get_max_prefix_len());
		string s = r.s.substr(1);
		int L = s.size();
		int starts = -1;
		double pushed_rate = (double)aligner->read_counters.pushed.get() / L;
		double popped_rate = (double)aligner->read_counters.popped.get() / L;
		double repeat_rate = (double)aligner->_repeated_visits / aligner->read_counters.pushed.get();
		double astar_missrate = 100.0 * astar.get_cache_misses() / astar.get_cache_trees();
		*pushed_rate_sum += pushed_rate;
		*popped_rate_sum += popped_rate;
		*repeat_rate_sum += repeat_rate;
		*pushed_rate_max = max(*pushed_rate_max, pushed_rate);
		*popped_rate_max = max(*popped_rate_max, popped_rate);
		*repeat_rate_max = max(*repeat_rate_max, repeat_rate);

		line[0] = 0;
		sprintf(line,
				"%8s\t%3d\t"
				"%8s\t%3d\t%d\t%4lf\t"
				"%8s\t%15s\t%8lf\t"
				"%3d\t%10s\t%10s\t"
				"%d\t%6d\t%6lf\t"
				"%6lf\t%4lf\t%8lf\t"
				"%8lf\t%2lf\t%d\n",
				args.graph_file, (int)aligner->graph().nodes(),
				algo.c_str(), astar.get_max_prefix_len(), int(astar.get_max_prefix_cost()), 100.0 * astar.get_compressable_vertices() / aligner->graph().nodes(),
				precomp_str.c_str(), r.comment.c_str(), 0.0,
				L, s.c_str(), spell(*path).c_str(),
				int(ans.cost), starts, pushed_rate,
				popped_rate, repeat_rate, timers.total.get_sec(),
				timers.astar.get_sec(), astar_missrate, aligner->unique_best);
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
	cout << "----" << endl;

	Measurers T;
	dict_t stats;   // string key -> string value

	LOG_DEBUG << "memory " << MemoryMeasurer::get_mem_gb();
	T.total.start();
    auto start_wt = std::chrono::high_resolution_clock::now();

	args = read_args(argc, argv);
	std::ios_base::sync_with_stdio(false);

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

	LOG_INFO << " ------------------------------ ";
	LOG_INFO << "Starting " << to_str(argc, argv);
	LOG_INFO << " ------------------------------ ";
	LOG_INFO << " sizeof(edge_t) = " << sizeof(edge_t);

	if (!performance_file.empty()) {
		FILE *fout = fopen(performance_file.c_str(), "w");
		fprintf(fout, "ref\trefsize\t"
				"algo\tA*-len-cap\tA*-cost-cap\tA*-compressable-vertices\t"
				"operation\treadname\tmemory\t"
				"len\tread\tspell\t"
				"cost\tstarts\tpushed\t"
				"popped\trepeat_rate\tt(map)\t"
				"t(astar)\tastar_missrate\tunique_best\n");
		fclose(fout);
	}

    graph_t G;
	vector<read_t> R;
	clock_t start;

	cout << "Loading reference graph... " << flush;
	T.read_graph.start();
	read_graph(&G, args.graph_file, output_dir);
	T.read_graph.stop();
	cout << "done." << endl << flush;

	cout << "Loading queries... " << flush;
	T.read_queries.start();
	read_queries(args.query_file, &R);
	T.read_queries.stop();
	cout << "done." << endl << flush;

	auto_params(G, R, &args);

	cout << "Contructing trie... " << flush;
	T.construct_trie.start();
	add_tree(&G, args.tree_depth);
	T.construct_trie.stop();
	cout << "done." << endl << flush;

	T.precompute.start();
	AStar astar(G, args.costs, args.AStarLengthCap, args.AStarCostCap, args.AStarNodeEqivClasses);
	T.precompute.stop();

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
		out << "             Alignment algo: " << args.algorithm 			 							<< endl;
		out << "                 Edit costs: " << int(args.costs.match) << ", " << int(args.costs.subst) << ", "
									<< int(args.costs.ins) << ", " << int(args.costs.del) << " (match, subst, ins, del)" << endl;
		out << "              Greedy match?: " << bool2str(args.greedy_match) 							<< endl;
		out << "                    Threads: " << args.threads											<< endl;
		out << " == A* parameters =="																<< endl;
		out << "                   Cost cap: " << args.AStarCostCap 									<< endl;
		out << "   Upcoming seq. length cap: " << args.AStarLengthCap 									<< endl;
		out << "      Nodes equiv. classes?: " << bool2str(args.AStarNodeEqivClasses) 					<< endl;
		out << "A* compressible equiv nodes: " << astar.get_compressable_vertices()
							<< " (" << 100.0 * astar.get_compressable_vertices() / G.nodes() << "%)"	<< endl;

		out << " == Data =="                                                                      << endl;
		// Note: the trie is built on top of the **doubled** original graph (incl. reverse).
		out << "         Original reference: " << G.orig_nodes << " nodes, " << G.orig_edges << " edges"<< endl;
		out << "                       Trie: " << G.trie_nodes << " nodes, " << G.trie_edges << " edges, "
												<< "depth: " << args.tree_depth 						<< endl;
		out << "  Reference+ReverseRef+Trie: " << G.nodes() << " nodes, " << G.edges() << " edges, "
												<< "density: " << (G.edges() / 2) / (G.nodes() / 2 * G.nodes() / 2) << endl;
		out << "                      Reads: " << R.size() << " x " << R.front().s.size()-1 << "bp, "
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
	std::mutex timer_m;
	Counter popped_trie_total, popped_ref_total;
	std::atomic<bool> nonunique_best_alignments(0);

	cout << "Aligning..." << flush;
	bool calc_mapping_cost = false;
	if (args.threads == 1) {
		FILE *fout = fopen(performance_file.c_str(), "a");
		for (size_t i=0; i<R.size(); i++) {
			char line[10000];
			Aligner aligner(G, align_params, &astar);
			state_t a_star_ans = wrap_readmap(R[i], algo, performance_file, &aligner, calc_mapping_cost,
					&R[i].edge_path, &pushed_rate_sum, &popped_rate_sum, &repeat_rate_sum, &pushed_rate_max, &popped_rate_max, &repeat_rate_max, line);
			fprintf(fout, "%s", line);
			total_timers += aligner.read_timers;
			if (!aligner.unique_best)
				nonunique_best_alignments = nonunique_best_alignments+1;

			if (i % (R.size() / 10) == 0) {
				cout << "A*-memoization at " << 100.0 * i / R.size() << "% of the reads aligned"
				<< ", entries: " << astar.entries() << ", "
				<< 100.0*b2gb(astar.table_mem_bytes_lower()) / MemoryMeasurer::get_mem_gb() << "%-"
				<< 100.0*b2gb(astar.table_mem_bytes_upper()) / MemoryMeasurer::get_mem_gb() << "%" << endl;
			}
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
				LOG_INFO << "thread " << t << " for reads [" << from << ", " << to << ")";
				for (size_t i=from; i<to; i++) {
					char line[10000];
					Aligner aligner(G, align_params, &astar);
					state_t a_star_ans = wrap_readmap(R[i], algo, performance_file, &aligner, calc_mapping_cost,
							&R[i].edge_path, &pushed_rate_sum, &popped_rate_sum, &repeat_rate_sum, &pushed_rate_max, &popped_rate_max, &repeat_rate_max, line);
					profileQueue.enqueue(string(line));
					{
						timer_m.lock();
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
	cout << "done." << endl << flush;

	{
		double astar_missrate = 100.0 * astar.get_cache_misses() / astar.get_cache_trees();
		double total_map_time = total_timers.total.get_sec();
		double total_mem = MemoryMeasurer::get_mem_gb();

		out << " == Aligning statistics =="														<< endl;
		out << "   Explored rate (avg, max): " << pushed_rate_sum / R.size() << ", " << pushed_rate_max << "    [states/bp] (states normalized by query length)" << endl;
		out << "     Expand rate (avg, max): " << popped_rate_sum / R.size() << ", " << popped_rate_max << endl;
		out << "      Memoization miss rate: " << astar_missrate << "%"									<< endl;
		out << "             Average popped: " << 1.0 * popped_trie_total.get() / (R.size()/args.threads) << " from trie vs "
											<< 1.0 * popped_ref_total.get() / (R.size()/args.threads) << " from ref"
											<< " (influenced by greedy match flag)" << endl;
		out << " Non-unique best alignments: " << nonunique_best_alignments << " out of " << R.size()	<< endl;

		out << " == Performance =="																<< endl;
		out << "    Memory: " << "                   measured | estimated" 									<< endl;
		out << "                   total: " << total_mem << "gb, 100% | -"		<< endl;
		out << "               reference: " << T.read_graph.m.get_gb() << "gb, " << 100.0*T.read_graph.m.get_gb() / total_mem << "% | " << 100.0*b2gb(G.reference_mem_bytes()) / total_mem << "%" << endl;
		out << "                   reads: " << T.read_queries.m.get_gb() << "gb, " << 100.0*T.read_queries.m.get_gb() / total_mem << "% | " << 100.0*b2gb(R.size() * R.front().size()) / total_mem << "%" << endl;
		out << "                    trie: " << T.construct_trie.m.get_gb() << "gb, " << 100.0*T.construct_trie.m.get_gb() / total_mem << "% | " << 100.0*b2gb(G.trie_mem_bytes()) / total_mem << "%" << endl;
		out << "     equiv. classes opt.: " << T.precompute.m.get_gb() << "gb, " << 100.0*T.precompute.m.get_gb() / total_mem << "% | " << 100.0*b2gb(astar.equiv_classes_mem_bytes()) / total_mem << "%" << endl;
		out << "          A*-memoization: " << T.align.m.get_gb() << "gb, " << 100.0*T.align.m.get_gb() / total_mem << "% | "
				     			<< 100.0*b2gb(astar.table_mem_bytes_lower()) / total_mem << "%-" 
				     			<< 100.0*b2gb(astar.table_mem_bytes_upper()) / total_mem << "%"
				     			<< " (" << int(astar.table_entrees()) << " entries)" 	<< endl;
		out << "    Wall runtime: " << total_wt.count() << "s"							<< endl;
		out << "       reference loading: " << T.read_graph.t.get_sec() << "s" 			<< endl;
		out << "         queries loading: " << T.read_queries.t.get_sec() << "s"			<< endl;
		out << "          construct trie: " << T.construct_trie.t.get_sec() << "s"		<< endl;
		out << "              precompute: " << T.precompute.t.get_sec() << "s"			<< endl;
		out << "                   align: " << align_wt.count() << "s (wall time) => "
				     						<< R.size() / total_map_time << " reads/s <=> "
				     						<< size_sum(R) / 1024.0 / total_map_time << " Kbp/s"	<< endl; 
		out << "                          " << T.align.t.get_sec() << "s (cpu time)"
				     << " (A*: " << 100.0 * total_timers.astar.get_sec() / total_map_time	<< "%, "
				     << "que: " << 100.0 * total_timers.queue.get_sec() / total_map_time	<< "%, "
				     << "dicts: " << 100.0 * total_timers.dicts.get_sec() / total_map_time	<< "%, "
				     << "greedy_match: " << 100.0 * total_timers.ff.get_sec() / total_map_time	<< "%"
				     << ")" << endl;
		out << " DONE" << endl;
	}

	extract_args_to_dict(args, &stats);
	T.extract_to_dict(&stats);
	ofstream tsv(stats_file);
	print_tsv(stats, tsv);
	tsv.close();

    return 0;
}

