#pragma once

#include <algorithm>
#include <map>
#include <unordered_set>
#include <string>
#include <vector>

#include "graph.h"
#include "utils.h"
#include "io.h"

typedef int seed_t;

namespace astarix {

class AStarSeedsWithErrors: public AStarHeuristic {
    // A*-seeds params
  public:
    struct Args {
        int seed_len;

		bool linear_reference;			// assume that the reference is linear
		bool match_pos_optimization;	// use the position of seed matches in order to maximize the heuristic function
		bool skip_near_crumbs;			// put crumbs in the trie only for nodes within [-m-delta, -m+delta] instead of [-m-delta, 0]

		// NOT USED
        int max_seed_errors;
        enum backwards_algo_t { DFS_FOR_LINEAR, BFS, COMPLEX, TOPSORT };
        backwards_algo_t backwards_algo;
		bool interval_intersection;  // only counting crumbs with intersecting intervals
    };

  private:
    struct Stats {
        Counter<> seeds;                    // number of seeds (depends only on the read)
        Counter<> seed_matches;             // places in the graph where seeds match well enough
        Counter<> paths_considered;         // number of paths updated for all seeds (supersource --> match)
        Counter<> states_with_crumbs;       // the number of states with crumbs
        Counter<> repeated_states;          // number of times crumbs are put on a states that already has crumbs
        Counter<cost_t> root_heuristic;     // heuristic from the trie root
        Counter<> heuristic_potential;      // maximal possible heuristic
        Counter<> reads;                    // reads processed

        void clear() {
            seeds.clear();
            seed_matches.clear();
            paths_considered.clear();
            states_with_crumbs.clear();
            repeated_states.clear();
            root_heuristic.clear();
        }

        Stats& operator+=(const Stats &b) {
            seeds += b.seeds;
            seed_matches += b.seed_matches;
            paths_considered += b.paths_considered;
            states_with_crumbs += b.states_with_crumbs;
            repeated_states += b.repeated_states;
            root_heuristic += b.root_heuristic;
            heuristic_potential += b.heuristic_potential;
            reads += b.reads;
            return *this;
        }
    };
	
	struct crumb_t {
		node_t v;  // can be in the trie or not
		node_t match_v;  // cannot be in the trie
		seed_t s;  // 0 for the rightmost seed
		crumb_t(const node_t _v, const node_t _match_v, const seed_t _s)
			: v(_v), match_v(_match_v), s(_s) {}
		bool operator<(const crumb_t &other) const {
			if (v != other.v)
				return v < other.v;
			if (match_v != other.match_v)
				return match_v < other.match_v;
			return s < other.s;
		}
	};

    // Fixed parameters
    const graph_t &G;
    const read_t *r_;
    const EditCosts &costs;
    const Args args;

	// Read alignment state
	int seeds_;
    int max_indels_;
    //std::unordered_map< node_t, std::vector<std::pair<node_t, seed_t> > > C;                        // node -> { match_position -> seed } }
	std::vector<crumb_t> C;   // all crumbs sorted by node, then by matching node, then by seed

	//static constexpr double learning_rate_ = 0.1;
	double retain_frac_;

	// Stats
    Stats read_cnt, global_cnt;

  public:
    AStarSeedsWithErrors(const graph_t &_G, const EditCosts &_costs, const Args &_args)
        : G(_G), costs(_costs), args(_args) {
		retain_frac_ = 0.1;  // start with using all seeds
    }

    // Cut r into chunks of length seed_len, starting from the end.
    void before_every_alignment(const read_t *r) {
		assert(C.empty());
        r_ = r;  // potentially useful for after_every_alignment()
		LOG_DEBUG << "max_indels: " << max_indels_;

        read_cnt.clear();
        read_cnt.reads.set(1);

		std::vector<pos_t> seed_starts = generate_seeds(r, retain_frac_);
		seeds_ = seed_starts.size();    //seeds_ = r->len / args.seed_len;  // TODO: not all seeds
		max_indels_ = (r->len * costs.match + seeds_ * costs.get_delta_min_special()) / costs.del;  // fuller

		generate_seeds_match_put_crumbs(seed_starts, r);
		std::sort(C.begin(), C.end());

		read_cnt.states_with_crumbs.set(C.size());
        read_cnt.seeds.set(seeds_);
        read_cnt.root_heuristic.set( h(state_t(0.0, 0, 0, -1, -1)) );
        read_cnt.heuristic_potential.set(seeds_);
        log_read_stats();

        global_cnt += read_cnt;
    }

	cost_t h(const state_t &st) const {
		int seeds_to_end = std::min((r_->len - st.i - 1) / args.seed_len, seeds_);  // TODO: not all seeds
		int missing = seeds_to_end;  // the maximum number of errors

		const auto from = std::lower_bound(C.begin(), C.end(), crumb_t(st.v, 0, 0));
		const auto to = std::lower_bound(C.begin(), C.end(), crumb_t(st.v+1, 0, 0));

		std::vector<int> cnt(seeds_to_end);
		int m = missing;  // m is the current number of different seeds inside of [it1, it2]

		for (auto it1=from, it2=from; it1 != to && it2 != to; ++it1) {  // left iterator
			for (; it2 != to && it2->match_v - it1->match_v <= max_indels_ + r_->len; ++it2)  // right itarator
				if (it2->s < seeds_to_end)
					if (cnt[it2->s]++ == 0)
						missing = std::min(missing, --m);
			if (it1->s < seeds_to_end)
				if (--cnt[it1->s] == 0)
					++m;
		}

		//for (auto it1 = from; it1 != to; ++it1) {
		//	std::unordered_set<seed_t> S;
		//	for (auto it2 = it1; it2 != to; ++it2) {
		//		assert(it2->match_v >= it1->match_v);
		//		if (it2->match_v - it1->match_v > max_indels_ + r_->len)
		//			break;
		//		if (it2->s < seeds_to_end)
		//			S.insert(it2->s);
		//	}
		//	missing = std::min(missing, seeds_to_end - (int)S.size());
		//}

		//const auto crumbs_it = C.find(st.v);
		//if (crumbs_it != C.end()) {
		//	auto crumbs = std::map(crumbs_it->second.begin(), crumbs_it->second.end());  // convert vector<pair> to map
		//	for (auto it1 = crumbs.begin(); it1 != crumbs.end(); ++it1) {  // *it1 == [match_v, seed]: 
		//		std::unordered_set<seed_t> S;
		//		for (auto it2 = it1; it2 != crumbs.end(); ++it2) {
		//			assert(it2->first >= it1->first);
		//			if (it2->first - it1->first > max_indels_ + r_->len) // + BINNING)
		//				break;
		//			if (it2->second < seeds_to_end)
		//				S.insert(it2->second);
		//		}
		//		missing = std::min(missing, seeds_to_end - (int)S.size());
		//	}
		//}
   
		return (r_->len - st.i)*costs.match + missing*costs.get_delta_min_special(); // fuller
	}

	inline int sign(double x) { return x < 0.0 ? -1 : +1; }

    void after_every_alignment(const AlignerTimers &t) {
        C.clear();  // clean up all crumbs for the next alignment

		// modify the retain_frac to keep the precomputation+A*query to take 50% of the time
		//double astar_frac_time = (t.astar_prepare_reads.get_sec() + t.astar.get_sec()) / t.total.get_sec();
		//retain_frac_ += 0.3*retain_frac_*sign(0.5-astar_frac_time);    // update the number of seeds
		//if (retain_frac_ > 1.0) retain_frac_ = 1.0;
		//else if (retain_frac_ < 0.0) retain_frac_ = 0.0;

		//LOG_DEBUG << "astar_frac_time: " << astar_frac_time << ", retain_frac_: " << retain_frac_;
    }

  private:
    inline void add_crumb_to_node(const seed_t p, const node_t match_v, const node_t curr_v) {
		if (args.match_pos_optimization)
			C.push_back(crumb_t(curr_v,match_v,p));
		else
			C.push_back(crumb_t(curr_v,0,p));  // turned off match position optimization

		//auto crumbs_it = C.find(curr_v);
		//node_t end_v = match_v;  // + p*args.seed_len;
		//if (crumbs_it == C.end()) {
		//	++read_cnt.states_with_crumbs;
		//	C[curr_v] = std::vector< std::pair<node_t, seed_t> >{ {end_v, p} };
		//}
		//else {
		//	auto &seeds = crumbs_it->second;
		//	//if (seeds.contains(end_v)) {
		//	//	++read_cnt.repeated_states;    
		//	//} else {
		//		seeds.push_back(std::make_pair(end_v, p));
		//		++read_cnt.states_with_crumbs;
		//	//}
		//}
		////assert(C.contains(curr_v) && C[curr_v].contains(end_v));
    }

	// NOT USED. BFS solution
    void add_crumbs_backwards_bfs(const seed_t p, const node_t match_v, int i, const node_t curr_v) {
		std::queue<node_t> Q({curr_v});
		std::unordered_map<node_t, int> pos = { {curr_v, i} };  // from curr_v

		while (!Q.empty()) {
			node_t v = Q.front(); Q.pop();
			int i = pos[v];
			add_crumb_to_node(p, match_v, v);
			for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it) {
				if (!args.skip_near_crumbs || (!G.node_in_trie(it->to) || i-1 <= G.get_trie_depth() + max_indels_))  // +/-delta optimization; assuming trie_depth <= G.get_trie_depth()  // TODO: +max_indels is not needed; change to max_deletions; 
					if (i > -max_indels_)
						if (!pos.contains(it->to)) {
							pos[it->to] = i-1;
							Q.push(it->to);
						}
			}
		}
	}

	// NOT USED
	void add_crumbs_backwards_anothersolution() {
		// solution 1
		//int curr_arr=0;
		//std::unordered_set<node_t> S[2];
		//S[curr_arr].insert(curr_v);
		//for (; i > -max_indels_; --i, curr_arr^=1) {
		//	S[curr_arr^1].clear();
		//	for (const node_t v: S[curr_arr])
		//		for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it)
		//			S[curr_arr^1].insert(it->to);
		//}
		//for (const node_t v: S[curr_arr])
		//	add_crumb_to_node(p, match_v, v);
	}

    void update_crumbs_up_the_trie(const seed_t p, const node_t match_v, node_t v) {
		if (v != 0)
			for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it)
				if (G.node_in_trie(it->to)) {
					add_crumb_to_node(p, match_v, it->to);
					update_crumbs_up_the_trie(p, match_v, it->to);
				}
    }

	// Breath-First-TopSort
    void add_crumbs_backwards(const seed_t p, const node_t match_v, int i, const node_t curr_v) {
		std::unordered_map<node_t, int> min_dist;
		std::unordered_map<node_t, int> max_dist;
		std::unordered_map<node_t, int> outgoing;
		std::queue<node_t> Q;
		edge_t e;

		Q.push(match_v);
		min_dist[match_v] = i;
		max_dist[match_v] = i;
		while(!Q.empty()) {
			node_t v = Q.front(); Q.pop();
			assert(min_dist.contains(v));
			assert(max_dist.contains(v));
			add_crumb_to_node(p, match_v, v);
			if (min_dist[v] <= G.get_trie_depth() + max_indels_) {
				update_crumbs_up_the_trie(p, match_v, v);
			}
			for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it) {
				node_t u = it->to;
				if (!G.node_in_trie(u)) {
					if (outgoing.contains(u)) outgoing[u] = outgoing[u]+1;
					else {
						outgoing[u] = 1;
						assert(!max_dist.contains(u));
						max_dist[u] = max_dist[v]-1;
					}
					assert(outgoing.contains(u));
					if (G.numOutOrigEdges(u,&e) == outgoing[u]) {
						assert(max_dist.contains(u));
						if (max_dist[u] >= -max_indels_) {
							Q.push(u);
							assert(!min_dist.contains(u));
							min_dist[u] = min_dist[v] - 1;
							//outgoing.remove(u);
						}
					}
				}
			}
		}
		//for (const auto it: outgoing)
		//	Q.push();
	}

	// For linear case only
    // `C[u]++` for all nodes (incl. `v`) that lead from `0` to `v` with a path of length `i`
    // Fully ignores labels.
    // Returns if the the supersource was reached at least once.
	// TODO: match back the seed letters (make sure it is not exponential: maybe reimplement with sets as in the paper)
    void add_crumbs_backwards_linear(const seed_t p, const node_t match_v, int i, const node_t curr_v) {
		// solution 0: fastest
		add_crumb_to_node(p, match_v, curr_v);
		if (i > -max_indels_)
			for (auto it=G.begin_orig_rev_edges(curr_v); it!=G.end_orig_rev_edges(); ++it)
				if (!args.skip_near_crumbs || !G.node_in_trie(it->to) || i-1 <= G.get_trie_depth())  // +/-delta optimization
					add_crumbs_backwards(p, match_v, i-1, it->to);
    }

    // Assumes that seed_len <= D so there are no duplicating outgoing labels.
	// Aligns a seed to the trie and stops when it goes to the reference
    void match_seed_put_crumbs(const read_t *r, const seed_t p, const int start, const int i, const node_t v) {
        if (i < start + args.seed_len) {
            // Match exactly down the trie and then through the original graph.
            for (auto it=G.begin_all_matching_edges(v, r->s[i]); it!=G.end_all_edges(); ++it) {
                if (it->type == ORIG || it->type == JUMP)  // ORIG in the graph, JUMP in the trie
                    match_seed_put_crumbs(r, p, start, i+1, it->to);
            }
        } else {
			add_crumbs_backwards(p, v, i, v);   // DEBUG: first go back to the starting of the match
            ++read_cnt.seed_matches;  // debug info
        }
    }

    // Split r into seeds of length seed_len. Return a vector of starting points of the seeds from last to first;
	std::vector<pos_t> generate_seeds(const read_t *r, double retain_frac) {
		int max_seeds = (r->len + args.seed_len - 1) / args.seed_len;  // ceil of r->len/seed_len
		std::vector<pos_t> seed_starts;
        for (int i=r->len-args.seed_len; i>=0; i-=args.seed_len) {
			seed_starts.push_back(i);
			if ((int)seed_starts.size() == int(retain_frac * max_seeds))
				break;
		}
		return seed_starts;
	}

	
    // For each exact occurence of a seed (i,v) in the graph,
    //   add 1 to C[u] for all nodes u on the path of match-length exactly `i` from supersource `0` to `v`
    void generate_seeds_match_put_crumbs(const std::vector<pos_t> &seed_starts, const read_t *r) {
		for (int i=0; i<(int)seed_starts.size(); i++) {
			pos_t st = seed_starts[i];
            match_seed_put_crumbs(r, i, st, st, 0);  // seed from [st, st+seed_len)
		}
    }

  public:
    void print_params(std::ostream &out) const {
        out << "          seed length: " << args.seed_len << " bp"         << std::endl;
    }

    void log_read_stats() {
        LOG_INFO << r_->comment << " A* seeds stats: "
            << read_cnt.seeds.get() << " seeds " 
            << "matching at " << read_cnt.seed_matches << " graph positions "
            << "and generating " << read_cnt.paths_considered << " paths "
            << "over " << read_cnt.states_with_crumbs << " states"
            << "(" << read_cnt.repeated_states << " repeated)"
            << "with best heuristic " << read_cnt.root_heuristic.get() << " "
			<< "with " << max_indels_ << "max_indels"; }

    void print_stats(std::ostream &out) const {
        int reads = global_cnt.reads.get();

        out << "        For all reads:"                                                     << std::endl;
        out << "                            Seeds: " << global_cnt.seeds << " (" << 1.0*global_cnt.seeds.get()/reads << " per read)"                  << std::endl;
        out << "                     Seed matches: " << global_cnt.seed_matches << " (" << 1.0*global_cnt.seed_matches.get()/reads << " per read, " << 1.0*global_cnt.seed_matches.get()/global_cnt.seeds.get() << " per seed)" << std::endl;
		out << "           Retained seed fraction: " << retain_frac_ << std::endl;
        out << "                 Paths considered: " << global_cnt.paths_considered << " (" << 1.0*global_cnt.paths_considered.get()/reads << " per read)"     << std::endl;
        out << "               States with crumbs: " << global_cnt.states_with_crumbs
            << " [+" << 100.0*global_cnt.repeated_states.get()/(global_cnt.states_with_crumbs.get()+global_cnt.repeated_states.get()) << "% repeated], (" << 1.0*global_cnt.states_with_crumbs.get()/reads << " per read)" << std::endl;
        out << "                  Heuristic (avg): " << 1.0*global_cnt.root_heuristic.get()/reads << " of potential " << 1.0*global_cnt.heuristic_potential.get()/reads << std::endl;
    }

};

}
