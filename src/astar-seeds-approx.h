#pragma once

#include <algorithm>
#include <map>
#include <unordered_set>
#include <string>
#include <vector>

#include "graph.h"
#include "utils.h"
#include "io.h"

typedef int node_t;
typedef int seed_t;

namespace astarix {

class AStarSeedsWithErrors: public AStarHeuristic {
    // A*-seeds params
  public:
    struct Args {
        enum backwards_algo_t { DFS_FOR_LINEAR, BFS, COMPLEX, TOPSORT };

        int seed_len;
        int max_seed_errors;
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
		node_t v;
		node_t match_v;
		seed_t s;
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

	//const int BINNING = 10000;

	// Stats
    Stats read_cnt, global_cnt;

  public:
    AStarSeedsWithErrors(const graph_t &_G, const EditCosts &_costs, const Args &_args)
        : G(_G), costs(_costs), args(_args) {
    }

    // Cut r into chunks of length seed_len, starting from the end.
    void before_every_alignment(const read_t *r) {
		assert(C.empty());
        r_ = r;  // potentially useful for after_every_alignment()
		seeds_ = r->len / args.seed_len;
		max_indels_ = (r->len * costs.match + seeds_ * costs.get_delta_min_special()) / costs.del;  // fuller
		LOG_DEBUG << "max_indels: " << max_indels_;

        read_cnt.clear();
        read_cnt.reads.set(1);

		generate_seeds_match_put_crumbs(r);
		std::sort(C.begin(), C.end());

        read_cnt.seeds.set(seeds_);
        read_cnt.root_heuristic.set( h(state_t(0.0, 0, 0, -1, -1)) );
        read_cnt.heuristic_potential.set(seeds_);
        log_read_stats();

        global_cnt += read_cnt;
    }

	cost_t h(const state_t &st) const {
		int seeds_to_end = (r_->len - st.i - 1) / args.seed_len;
		int missing = seeds_to_end;  // the maximum number of errors

		const auto from = std::lower_bound(C.begin(), C.end(), crumb_t(st.v, 0, 0));
		const auto to = std::lower_bound(C.begin(), C.end(), crumb_t(st.v+1, 0, 0));

		std::vector<int> cnt(seeds_to_end);
		int m = missing;

		for (auto it1=from, it2=from; it1 != to && it2 != to; ++it1) {
			for (; it2 != to && it2->match_v - it1->match_v <= max_indels_ + r_->len; ++it2)
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

    void after_every_alignment() {
        C.clear();  // clean up all crumbs for the next alignment
    }

  private:
    inline void add_crumb_to_node(const seed_t p, const node_t match_v, const node_t curr_v) {
		C.push_back(crumb_t(curr_v,match_v,p));

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

    // `C[u]++` for all nodes (incl. `v`) that lead from `0` to `v` with a path of length `i`
    // Fully ignores labels.
    // Returns if the the supersource was reached at least once.
	// TODO: match back the seed letters (make sure it is not exponential: maybe reimplement with sets as in the paper)
    void add_crumbs_backwards(const seed_t p, const node_t match_v, int i, const node_t curr_v) {
		// solution 2
		std::queue<node_t> Q({curr_v});
		std::unordered_map<node_t, int> pos = { {curr_v, i} };  // from curr_v

		//while (!Q.empty()) {
		//	node_t v = Q.front(); Q.pop();
		//	add_crumb_to_node(p, match_v, v);
		//	if (pos[v] <= -max_indels_)
		//		break;
		//	for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it) {
		//		//if (G.node_in_trie(it->to) && pos[v] > args.seed_len + max_indels_)  // assuming trie_depth <= seed_len
		//		//	continue;
		//		if (!pos.contains(it->to)) {
		//			pos[it->to] = pos[v]-1;
		//			Q.push(it->to);
		//		}
		//	}
		//}

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

		// solution 0: fastest
		add_crumb_to_node(p, match_v, curr_v);
		if (i > -max_indels_)
			for (auto it=G.begin_orig_rev_edges(curr_v); it!=G.end_orig_rev_edges(); ++it)
				if (!G.node_in_trie(it->to) || i-1 <= args.seed_len + max_indels_)  // assuming trie_depth <= seed_len
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
            assert(!G.node_in_trie(v));
			add_crumbs_backwards(p, v, i, v);   // DEBUG: first go back to the starting of the match
            ++read_cnt.seed_matches;  // debug info
        }
    }

    // Split r into seeds of length seed_len.
    // For each exact occurence of a seed (i,v) in the graph,
    //   add 1 to C[u] for all nodes u on the path of match-length exactly `i` from supersource `0` to `v`
    // Returns the number of seeds.
    int generate_seeds_match_put_crumbs(const read_t *r) {
        int seeds = 0;
        for (int i=r->len-args.seed_len; i>=0; i-=args.seed_len) {
            match_seed_put_crumbs(r, seeds, i, i, 0);  // seed from [i, i+seed_len)
            seeds++;
        }
        return seeds;
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
			<< "with " << max_indels_ << "max_indels";
    }

    void print_stats(std::ostream &out) const {
        int reads = global_cnt.reads.get();

        out << "        For all reads:"                                                     << std::endl;
        out << "                            Seeds: " << global_cnt.seeds << " (" << 1.0*global_cnt.seeds.get()/reads << " per read)"                  << std::endl;
        out << "                     Seed matches: " << global_cnt.seed_matches << " (" << 1.0*global_cnt.seed_matches.get()/reads << " per read, " << 1.0*global_cnt.seed_matches.get()/global_cnt.seeds.get() << " per seed)" << std::endl;
        out << "                 Paths considered: " << global_cnt.paths_considered << " (" << 1.0*global_cnt.paths_considered.get()/reads << " per read)"     << std::endl;
        out << "               States with crumbs: " << global_cnt.states_with_crumbs
            << " [+" << 100.0*global_cnt.repeated_states.get()/(global_cnt.states_with_crumbs.get()+global_cnt.repeated_states.get()) << "% repeated], (" << 1.0*global_cnt.states_with_crumbs.get()/reads << " per read)" << std::endl;
        out << "                  Heuristic (avg): " << 1.0*global_cnt.root_heuristic.get()/reads << " of potential " << 1.0*global_cnt.heuristic_potential.get()/reads << std::endl;
    }

};

}
