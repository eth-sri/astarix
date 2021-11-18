#pragma once

#include <algorithm>
#include <map>
#include <unordered_set>
#include <string>
#include <vector>

#include "graph.h"
#include "utils.h"
#include "io.h"

typedef short seed_t;

namespace astarix {

class AStarSeedsWithErrors: public AStarHeuristic {
  public:
    // A*-seeds parameters.
    struct Args {
        int seed_len;
		bool skip_near_crumbs;			// Put crumbs in the trie only for nodes within [-m-delta, -m+delta] instead of [-m-delta, 0].
    };

  private:
    struct Stats {
        Counter<> seeds;                    // number of seeds (depends only on the read)
        Counter<> seed_matches;             // places in the graph where seeds match well enough
        Counter<> states_with_crumbs;       // the number of states with crumbs
        Counter<> repeated_states;          // number of times crumbs are put on a states that already has crumbs
        Counter<cost_t> root_heuristic;     // heuristic from the trie root
        Counter<> heuristic_potential;      // maximal possible heuristic
        Counter<> reads;                    // reads processed

        void clear() {
            seeds.clear();
            seed_matches.clear();
            states_with_crumbs.clear();
            repeated_states.clear();
            root_heuristic.clear();
        }

        Stats& operator+=(const Stats &b) {
            seeds += b.seeds;
            seed_matches += b.seed_matches;
            states_with_crumbs += b.states_with_crumbs;
            repeated_states += b.repeated_states;
            root_heuristic += b.root_heuristic;
            heuristic_potential += b.heuristic_potential;
            reads += b.reads;
            return *this;
        }
    };
	
	struct crumb_t {
		node_t v;  // Can be both in the trie or not.
		seed_t s;  // From right to left: the last/rightmost seed has index 0.

		crumb_t() {}
		crumb_t(const node_t _v, const seed_t _s)
			: v(_v), s(_s) {}

		bool operator<(const crumb_t &other) const {
			if (v != other.v)
				return v < other.v;
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
	std::vector<crumb_t> C;   // All crumbs sorted by node, then by seed number.

	// Stats
    Stats read_cnt, global_cnt;

  private:
	// Skip v (which is in the reference) and add crumbs only to the trie.
    void update_crumbs_up_the_trie(const seed_t s, const node_t match_v, node_t v) {
		if (v != 0)
			for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it)
				if (G.node_in_trie(it->to)) {
					add_crumb_to_node(s, match_v, it->to);
					update_crumbs_up_the_trie(s, match_v, it->to);
				}
    }

    inline void add_crumb_to_node(const seed_t s, const node_t match_v, const node_t curr_v) {
		C.push_back(crumb_t(curr_v,s)); 
    }

  public:
    AStarSeedsWithErrors(const graph_t &_G, const EditCosts &_costs, const Args &_args)
        : G(_G), costs(_costs), args(_args) {
		if (args.seed_len == -1)
			throw "seed len not set.";
    }

    // Cut r into chunks of length seed_len, starting from the end.
    void before_every_alignment(const read_t *r) {
		assert(C.empty());
        r_ = r;

        read_cnt.clear();
        read_cnt.reads.set(1);

		std::vector<pos_t> seed_starts = generate_seeds(r, 1.0);
		seeds_ = seed_starts.size(); 
		max_indels_ = std::ceil((r->len * costs.match + seeds_ * costs.get_delta_min_special()) / costs.del);
		LOG_DEBUG << "max_indels: " << max_indels_;

		match_all_seeds(seed_starts, r);
		std::sort(C.begin(), C.end());

		read_cnt.states_with_crumbs.set(C.size());
        read_cnt.seeds.set(seeds_);
        read_cnt.root_heuristic.set( h(state_t(0.0, 0, 0, -1, -1)) );
        read_cnt.heuristic_potential.set(seeds_);
        log_read_stats();

        global_cnt += read_cnt;
    }

    // Split r into seeds of length seed_len. Return a vector of starting points of the seeds from last to first;
	std::vector<pos_t> generate_seeds(const read_t *r, double seeds_retain_frac) {
		int max_seeds = (r->len + args.seed_len - 1) / args.seed_len;  // ceil of r->len/seed_len
		std::vector<pos_t> seed_starts;
        for (int i=r->len-args.seed_len; i>=0; i-=args.seed_len) {
			seed_starts.push_back(i);
			if ((int)seed_starts.size() == int(seeds_retain_frac * max_seeds))
				break;
		}
		return seed_starts;
	}
	
    // For each exact occurence of a seed (i,v) in the graph,
    //   add 1 to C[u] for all nodes u on the path of match-length exactly `i` from supersource `0` to `v`
    void match_all_seeds(const std::vector<pos_t> &seed_starts, const read_t *r) {
		for (int s=0; s<(int)seed_starts.size(); s++) {
			pos_t start_pos = seed_starts[s];
            match_reverse_complement_seed(r, s, start_pos, start_pos+args.seed_len-1, G.trie_root());  // seed from [st, st+seed_len)
		}
    }
	
    // Assumes that seed_len >= D so the match is outside of the trie.
    void match_reverse_complement_seed(const read_t *r, const seed_t s, const int start, const int i, const node_t v) {
        if (i >= start) {
            // Match exactly down the trie and then through the original graph.
			label_t c = compl_nucl(r->s[i]);
            for (auto it=G.begin_all_matching_edges(v, c); it!=G.end_all_edges(); ++it)
                if (it->type == ORIG || it->type == JUMP)  // ORIG in the graph, JUMP in the trie
                    match_reverse_complement_seed(r, s, start, i-1, it->to);
        } else {
			// All the seed is aligned now.
			node_t u = G.node2revcompl(v);
			put_crumbs_backwards(s, u, i);
            ++read_cnt.seed_matches;
        }
    }

	// TopSort from match_v on backwards edges with max distance i+max_indels_.
    void put_crumbs_backwards(const seed_t s, const node_t match_v, int i) {
		std::unordered_map<node_t, int> min_pos;                        // _minimal_ read index where an _expanded_ node can be aligned without indels so that r[i] aligns at match_v
		std::unordered_map<node_t, int> max_pos;                        // _maximal_ read index where an _explored_ node --||--
		std::unordered_map<node_t, int> outgoing;                       // Number of explored outgoing edges of a node
		std::queue<node_t> Q;
		edge_t e;

		bool start_in_a_loop = false;									// Handle the case when match_v is in a (small enough) cycle.

		// TopSort in referece (w/o trie)
		Q.push(match_v);
		min_pos[match_v] = max_pos[match_v] = i;
		while(!Q.empty() && !start_in_a_loop) {
			node_t v = Q.front(); Q.pop();
																		assert(min_pos.contains(v));
																		assert(max_pos.contains(v));
			add_crumb_to_node(s, match_v, v);
			if (!args.skip_near_crumbs || min_pos[v] <= G.get_trie_depth() + max_indels_)
				update_crumbs_up_the_trie(s, match_v, v);
			for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it) {
				node_t u = it->to;
				if (!G.node_in_trie(u)) {
					if (u == match_v) { 								// Looping back to the starting node.
						start_in_a_loop = true;
						break;
					}
					if (outgoing.contains(u)) outgoing[u] = outgoing[u]+1;
					else {
						outgoing[u] = 1;
																		assert(!max_pos.contains(u));
						max_pos[u] = max_pos[v]-1;
					}
																		assert(outgoing.contains(u));
					if (G.numOutOrigEdges(u,&e) == outgoing[u]) {
																		assert(max_pos.contains(u));
																		assert(!min_pos.contains(u));
						min_pos[u] = min_pos[v] - 1;
						// max_pos[u] = match_start_pos - mindist(u,v)  =>  dist(u,v) < i+n_del <=> max_pos[u] > -n_del
						if (max_pos[u] > -max_indels_) {
							Q.push(u);
						}
					}
				}
			}
		}

		if (start_in_a_loop) {
			// Start from the beginning. Does not matter if any crumbs were already added.
			max_pos.clear();
			max_pos[match_v] = i;
			Q.push(match_v);
		} else {
			// Initialize Q with nodes from explored but not expanded nodes (aka from cycles).
			for (const auto &[v, _]: max_pos)
				if (!min_pos.contains(v))
					Q.push(v);
			LOG_DEBUG << "Cycles reached: " << Q.size();
		}

		// BFS on both reference graph and trie: add a crumb to all nodes before position -max_indels_.
		while (!Q.empty()) {
			node_t v = Q.front(); Q.pop();
			add_crumb_to_node(s, match_v, v);
			for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it) {
				node_t u = it->to;
				if (!max_pos.contains(u))
					if (max_pos[v]-1 >= -max_indels_) {
						max_pos[u] = max_pos[v]-1;
						Q.push(u);
					}
			}
		}
	}

	// Seed heuristic query called during A* alignment.
	cost_t h(const state_t &st) const {
		int seeds_to_end = std::min((r_->len - st.i - 1) / args.seed_len, seeds_);
		int missing = seeds_to_end;  // Maximum number of errors.

		const auto from = std::lower_bound(C.begin(), C.end(), crumb_t(st.v, 0));
		const auto to = std::lower_bound(C.begin(), C.end(), crumb_t(st.v+1, 0));

		std::vector<int> cnt(seeds_to_end);
		int m = missing;  // m is the current number of different seeds inside of [it1, it2]

		for (auto it1=from, it2=from; it1 != to && it2 != to; ++it1) {  // left iterator
			for (; it2 != to; ++it2)
				if (it2->s < seeds_to_end)
					if (cnt[it2->s]++ == 0)
						missing = std::min(missing, --m);
			if (it1->s < seeds_to_end)
				if (--cnt[it1->s] == 0)
					++m;
		}
   
		return (r_->len - st.i)*costs.match + missing*costs.get_delta_min_special();
	}

    void after_every_alignment(const AlignerTimers &t) {
        C.clear();  // Clean up all crumbs before next alignment.
    }

    void print_params(std::ostream &out) const {
        out << "          seed length: " << args.seed_len << " bp"         << std::endl;
        out << "     skip near crumbs: " << args.skip_near_crumbs          << std::endl;
    }

    void print_stats(std::ostream &out) const {
        int reads = global_cnt.reads.get();

        out << "        For all reads:"                                                     << std::endl;
        out << "                            Seeds: " << global_cnt.seeds << " (" << 1.0*global_cnt.seeds.get()/reads << " per read)"                  << std::endl;
        out << "                     Seed matches: " << global_cnt.seed_matches << " (" << 1.0*global_cnt.seed_matches.get()/reads << " per read, " << 1.0*global_cnt.seed_matches.get()/global_cnt.seeds.get() << " per seed)" << std::endl;
        out << "               States with crumbs: " << global_cnt.states_with_crumbs
            << " [+" << 100.0*global_cnt.repeated_states.get()/(global_cnt.states_with_crumbs.get()+global_cnt.repeated_states.get()) << "% repeated], (" << 1.0*global_cnt.states_with_crumbs.get()/reads << " per read)" << std::endl;
        out << "                  Heuristic (avg): " << 1.0*global_cnt.root_heuristic.get()/reads << " of potential " << 1.0*global_cnt.heuristic_potential.get()/reads << std::endl;
    }

    void log_read_stats() {
        LOG_INFO << r_->comment << " A* seeds stats: "
            << read_cnt.seeds.get() << " seeds " 
            << "matching at " << read_cnt.seed_matches << " graph positions "
            << "over " << read_cnt.states_with_crumbs << " states"
            << "(" << read_cnt.repeated_states << " repeated)"
            << "with best heuristic " << read_cnt.root_heuristic.get() << " "
			<< "with " << max_indels_ << "max_indels";
	}

	int crumbs() const {
		return global_cnt.states_with_crumbs.get();
	}
};

}
