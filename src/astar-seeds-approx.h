#pragma once

#include <map>
#include <unordered_set>
#include <string>
#include <vector>

#include "graph.h"
#include "utils.h"
#include "io.h"

typedef int node_t;

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

    // Fixed parameters
    const graph_t &G;
    const read_t *r_;
    const EditCosts &costs;
    const Args args;

	// Read alignment state
    int max_indels_;
    std::unordered_map< node_t, std::unordered_set<int> > C;                        // node -> seeds

	// Stats
    Stats read_cnt, global_cnt;

  public:
    AStarSeedsWithErrors(const graph_t &_G, const EditCosts &_costs, const Args &_args)
        : G(_G), costs(_costs), args(_args) {
    }

    // Cut r into chunks of length seed_len, starting from the end.
    void before_every_alignment(const read_t *r) {
        r_ = r;  // potentially useful for after_every_alignment()

		int seeds = r->len / args.seed_len;
		max_indels_ = (r->len * costs.match + seeds * costs.get_delta_min_special()) / costs.del;  // fuller
		LOG_DEBUG << "max_indels: " << max_indels_;

        read_cnt.clear();
        read_cnt.reads.set(1);

		generate_seeds_match_put_crumbs(r);

        read_cnt.seeds.set(seeds);
        read_cnt.root_heuristic.set( h(state_t(0.0, 0, 0, -1, -1)) );
        read_cnt.heuristic_potential.set(seeds);
        log_read_stats();

        global_cnt += read_cnt;
    }

	cost_t h(const state_t &st) const {
		int seeds_to_end = (r_->len - st.i - 1) / args.seed_len;
		int missing = seeds_to_end;  // the maximum number of errors
				  
		const auto crumbs_it = C.find(st.v);
		if (crumbs_it != C.end())
			for (const auto &seed: crumbs_it->second)
				if (seed < seeds_to_end)  // the seeds are numbered from 0 starting from the back of the read
					--missing;
   
		return (r_->len - st.i)*costs.match + missing*costs.get_delta_min_special(); // fuller
	}

    void after_every_alignment() {
        C.clear();  // clean up all crumbs for the next alignment
    }

  private:
    inline void add_crumb_to_node(int p, node_t v) {
		auto crumbs_it = C.find(v);
		if (crumbs_it == C.end()) {
			++read_cnt.states_with_crumbs;

			std::unordered_set<int> tmp;
			tmp.insert(p);
			C[v] = tmp;
		}
		else {
			auto &seeds = crumbs_it->second;
			if (seeds.contains(p)) {
				++read_cnt.repeated_states;    
			} else {
				seeds.insert(p);
				++read_cnt.states_with_crumbs;
			}
		}
		assert(C.contains(v) && C[v].contains(p));
    }

    // `C[u]++` for all nodes (incl. `v`) that lead from `0` to `v` with a path of length `i`
    // Fully ignores labels.
    // Returns if the the supersource was reached at least once.
	// TODO: match back the seed letters (make sure it is not exponential: maybe reimplement with sets as in the paper)
    bool add_crumbs_backwards(int p, int i, node_t v) {
		add_crumb_to_node(p, v);
		if (i > -max_indels_)
			for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it)
				add_crumbs_backwards(p, i-1, it->to);

        return true;
    }

    // Assumes that seed_len <= D so there are no duplicating outgoing labels.
	// Aligns a seed to the trie and stops when it goes to the reference
    void match_seed_put_crumbs(const read_t *r, int p, int start, int i, node_t v) {
        if (i < start + args.seed_len) {
            // Match exactly down the trie and then through the original graph.
            for (auto it=G.begin_all_matching_edges(v, r->s[i]); it!=G.end_all_edges(); ++it) {
                if (it->type == ORIG || it->type == JUMP)  // ORIG in the graph, JUMP in the trie
                    match_seed_put_crumbs(r, p, start, i+1, it->to);
            }
        } else {
            assert(!G.node_in_trie(v));
			add_crumbs_backwards(p, i, v);   // DEBUG: first go back to the starting of the match
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
