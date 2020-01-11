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
		edge_path_t *path, double memory, double *pushed_rate_sum, double *popped_rate_sum, double *repeat_rate_sum, double *pushed_rate_max, double *popped_rate_max, double *repeat_rate_max, string *perf_s, char *line) {
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
		double astar_hitrate = 100.0 - 100.0 * astar.get_cache_misses() / astar.get_cache_trees();
		*pushed_rate_sum += pushed_rate;
		*popped_rate_sum += popped_rate;
		*repeat_rate_sum += repeat_rate;
		*pushed_rate_max = max(*pushed_rate_max, pushed_rate);
		*popped_rate_max = max(*popped_rate_max, popped_rate);
		*repeat_rate_max = max(*repeat_rate_max, repeat_rate);

		//char line[10000];
		line[0] = 0;
		sprintf(line,
				"%8s\t%3d\t"
				"%8s\t%3d\t%4lf\t%4lf\t"
				"%8s\t%15s\t%8lf\t"
				"%3d\t%10s\t%10s\t"
				"%8lf\t%6d\t%6lf\t"
				"%6lf\t%4lf\t%8lf\t"
				"%8lf\t%2lf\n",
				args.graph_file, (int)aligner->graph().nodes(),
				algo.c_str(), astar.get_max_prefix_len(), astar.get_max_prefix_cost(), 100.0 * astar.get_compressable_vertices() / aligner->graph().nodes(),
				precomp_str.c_str(), r.comment.c_str(), memory,
				L, s.c_str(), spell(*path).c_str(),
				ans.cost, starts, pushed_rate,
				popped_rate, repeat_rate, timers.total.get_sec(),
				timers.astar.get_sec(), astar_hitrate);
		//*perf_s += line;
		
		//if (perf_s->length() > 50000) {
		//	FILE *fout = fopen(performance_file.c_str(), "a");
		//	fprintf(fout, "%s", perf_s->c_str());
		//	fclose(fout);

		//	perf_s->clear();
		//}
	}

	return ans;
}

void read_queries(const char *query_file, vector<read_t> *R) {
	cout << "Reading queries..." << endl << flush;

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
		args->tree_depth = ceil(log(G.nodes()) / log(4.0));
	}
}

void print_tsv(map<string, string> dict) {
	for (auto const &it: dict)
		cout << it.first << "\t";
	cout << endl;
	for (auto const &it: dict) {
		assert(!it.second.empty());
		cout << it.second << "\t";
	}
	cout << endl;
}

typedef map<string, string> dict_t;

struct Timers {
	Timer input, align, read_graph, read_queries, precompute;

	void append_to_dict(dict_t *dict) {
		(*dict)["input_sec"] = to_string(input.get_sec());
		(*dict)["align_sec"] = to_string(align.get_sec());
	}
};

Timers T;
dict_t stats;   // string key -> string value

int main(int argc, char **argv) {
	T.input.start();

	args = read_args(argc, argv);
	std::ios_base::sync_with_stdio(false);

	string performance_file, info_log_file;

	string output_dir = args.output_dir;
	if (!output_dir.empty()) {
		assure_dir_exists(args.output_dir);
		performance_file = output_dir + "/alignments.tsv";
		info_log_file = output_dir + "/info.log";
		//stats_log_file = output_dir + "/info.log";
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
				"t(astar)\tastar_hitrate\n");
		fclose(fout);
	}

    graph_t G;
	vector<read_t> R;
	clock_t start;

	T.read_graph.start();
	read_graph(&G, args.graph_file, output_dir);
	T.read_graph.stop();

	T.read_queries.start();
	read_queries(args.query_file, &R);
	T.read_queries.stop();

	auto_params(G, R, &args);

	add_tree(&G, args.tree_depth);

	T.precompute.start();
	AStar astar(G, args.costs, args.AStarLengthCap, args.AStarCostCap, args.AStarNodeEqivClasses);
	T.precompute.stop();

	AlignParams align_params(args.costs, args.greedy_match); //, args.tree_depth);
	string algo = string(args.algorithm);
	string perf_s;

	assert(G.has_supersource());
	LOG_INFO << "Mapping init with graph with n=" << G.V.size() << " and m=" << G.E.size();
	align_params.print();

    double pushed_rate_sum(0.0), pushed_rate_max(0.0);
    double popped_rate_sum(0.0), popped_rate_max(0.0);
    double repeat_rate_sum(0.0), repeat_rate_max(0.0);
	//cost_t total_read_cost(0.0), max_read_cost(0.0);

    double precomp_vm, precomp_rss;
    process_mem_usage(precomp_vm, precomp_rss);
	precomp_vm /= 1024.0 * 1024.0;  // to GB
	precomp_rss /= 1024.0 * 1024.0;  // to GB

	T.input.stop();
	T.align.start();

	cout << "Aligning..." << endl << flush;
	bool calc_mapping_cost = false;
	if (false && args.threads == 1) {
		FILE *fout = fopen(performance_file.c_str(), "a");
		for (size_t i=0; i<R.size(); i++) {
			char line[10000];
			Aligner aligner(G, align_params, &astar);
			state_t a_star_ans = wrap_readmap(R[i], algo, performance_file, &aligner, calc_mapping_cost,
					&R[i].edge_path, precomp_vm, &pushed_rate_sum, &popped_rate_sum, &repeat_rate_sum, &pushed_rate_max, &popped_rate_max, &repeat_rate_max, &perf_s, line);
			fprintf(fout, "%s", line);
		}
		fclose(fout);
	} else {
		moodycamel::ConcurrentQueue<string> profileQueue;
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
							&R[i].edge_path, precomp_vm, &pushed_rate_sum, &popped_rate_sum, &repeat_rate_sum, &pushed_rate_max, &popped_rate_max, &repeat_rate_max, &perf_s, line);
					profileQueue.enqueue(string(line));
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

	std::ostream &out = cout;
	double astar_hitrate = 100.0 - 100.0 * astar.get_cache_misses() / astar.get_cache_trees();

//	double total_map_time = aligner.total_timers.total.get_sec();

	T.append_to_dict(&stats);

	out.setf(ios::fixed, ios::floatfield);
	out.precision(2);
	out << endl;
	out << "       == Input =="                                                                         << endl;
	out << "              Queries/reads: " << R.size() << " read in " << T.read_queries.get_sec() << "s"<< endl;
	out << "Graph nodes, edges, density: " << G.nodes() << ", " << G.edges()
													<< ", " << (G.edges() / 2) / (G.nodes() / 2 * G.nodes() / 2) << endl;
	out << "                    Threads: " << args.threads											<< endl;
	out << "           Reading run time: " << T.read_graph.get_sec() << "s" 	<< endl;  // There are the independent forward and reverse graph
	//out << "                  Read cost: " << "avg=" << total_read_cost / R.size()
//											<< ", max=" << max_read_cost							<< endl;
	out << "       == General parameters and optimizations == "                                     << endl;
	out << "             Alignment algo: " << args.algorithm 			 							<< endl;
	out << "                 Edit costs: " << args.costs.match << ", " << args.costs.subst << ", " << args.costs.ins << ", " << args.costs.del
		<< " (match, subst, ins, del)" << endl;
	out << "                    Threads: " << "1"													<< endl;
	out << "              Greedy match?: " << bool2str(args.greedy_match) 							<< endl;
	out << "       == Aligning statistics =="														<< endl;
	out << "                      Reads: " << R.size() << " x " << R.front().s.size()-1 << "bp"	<< endl;
	out << "         Aligning wall time: " << T.align.get_sec() << "s"								<< endl;
//								<< " (A*: " << 100.0 * aligner.total_timers.astar.get_sec() / total_map_time	<< "%, "
//								<< "que: " << 100.0 * aligner.total_timers.queue.get_sec() / total_map_time	<< "%, "
//								<< "dicts: " << 100.0 * aligner.total_timers.dicts.get_sec() / total_map_time	<< "%, "
//								<< "greedy_match: " << 100.0 * aligner.total_timers.ff.get_sec() / total_map_time	<< "%"
//								<< ")" << endl;
//	out << "                Performance: " << R.size() / total_map_time << " reads/s, "
//											<< size_sum(R) / 1024.0 / total_map_time << " Kbp/s"	<< endl;
	out << "   Explored rate (avg, max): " << pushed_rate_sum / R.size() << ", " << pushed_rate_max << "    [states/bp] (states normalized by query length)" << endl;
	out << "     Expand rate (avg, max): " << popped_rate_sum / R.size() << ", " << popped_rate_max << endl;

	if (algo == "astar-prefix") {
	out << "       == A* parameters =="																<< endl;
	out << "                   Cost cap: " << args.AStarCostCap 									<< endl;
	out << "   Upcoming seq. length cap: " << args.AStarLengthCap 									<< endl;
	out << "      Nodes equiv. classes?: " << bool2str(args.AStarNodeEqivClasses) 					<< endl;
	out << "  Compressible equiv. nodes: " << astar.get_compressable_vertices()
						<< " (" << 100.0 * astar.get_compressable_vertices() / G.nodes() << "%)"	<< endl;
	out << "    Precomputation run time: " << T.precompute.get_sec() << "s"							<< endl;
	out << "                     Memory: " << precomp_vm << " GB "
					    				<< "(graph: " << 100.0*G.GB() / precomp_vm << "%, " 	
										<< "A*-memoization: "
//										<< 100.0*astar.table_gb_lower() / precomp_vm << "%-"
										<< "<" << 100.0*astar.table_gb_upper() / precomp_vm << "%"
										<< " for " << int(astar.table_entrees()) << " entries)" << endl;
	out << "       Memoization hit rate: " << astar_hitrate << "%"									<< endl;
	}

	if (!performance_file.empty()) {
		FILE *fout = fopen(performance_file.c_str(), "a");
		fprintf(fout, "%s", perf_s.c_str());
		fclose(fout);
	}

	//print_tsv(stats);

    return 0;
}

