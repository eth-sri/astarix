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
    std::unordered_map< node_t, std::unordered_map<int, cost_t> > C;                        // node -> (seed->cost)

	// Stats
    Stats read_cnt, global_cnt;

  public:
    AStarSeedsWithErrors(const graph_t &_G, const EditCosts &_costs, const Args &_args)
        : G(_G), costs(_costs), args(_args) {
    }

    cost_t h(const state_t &st) const {
        int all_seeds_to_end = (r_->len - st.i - 1) / args.seed_len;
        int total_errors = (args.max_seed_errors+1)*all_seeds_to_end;  // the maximum number of errors
        
        const auto crumbs_it = C.find(st.v);
        if (crumbs_it != C.end())
            for (const auto &[seed, cost]: crumbs_it->second)
                if (seed < all_seeds_to_end)
                    total_errors -= cost;

        cost_t res = total_errors * costs.get_min_mismatch_cost(); // tested
        //cost_t res = (r_->len - st.i) * costs.match + total_errors * costs.get_delta_min_special(); // fuller

        return res;
    }

    // Cut r into chunks of length seed_len, starting from the end.
    void before_every_alignment(const read_t *r) {
        r_ = r;  // potentially useful for after_every_alignment()
        C.clear();

		int seeds = r->len / args.seed_len;
		//max_indels_ = (seeds * costs.get_min_mismatch_cost()) / costs.del;  // stable 
		max_indels_ = (r->len * costs.match + seeds * costs.get_delta_min_special()) / costs.del;  // fuller
		LOG_DEBUG << "max_indels: " << max_indels_;


        read_cnt.clear();
        read_cnt.reads.set(1);

		generate_seeds_match_put_crumbs(r);

        read_cnt.seeds.set(seeds);
        read_cnt.root_heuristic.set( h(state_t(0.0, 0, 0, -1, -1)) );
        read_cnt.heuristic_potential.set( (args.max_seed_errors+1)*read_cnt.seeds.get() );
        log_read_stats();

        global_cnt += read_cnt;
        states_with_crumbs = read_cnt.states_with_crumbs;  // TODO: refactor
    }

    void after_every_alignment() {
    }

    void print_params(std::ostream &out) const {
        out << "          seed length: " << args.seed_len << " bp"         << std::endl;
        out << "      max seed errors: " << args.max_seed_errors           << std::endl;
    }

    void log_read_stats() {
        LOG_INFO << r_->comment << " A* seeds stats: "
            << read_cnt.seeds.get() << " seeds " 
            << "matching at " << read_cnt.seed_matches << " graph positions "
            << "and generating " << read_cnt.paths_considered << " paths "
            << "over " << read_cnt.states_with_crumbs << " states"
            << "(" << read_cnt.repeated_states << " repeated)"
            << "with best heuristic " << read_cnt.root_heuristic.get() << " "
            << "out of possible " << (args.max_seed_errors+1)*read_cnt.seeds.get() << " "
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

  private:
    inline bool crumbs_already_set(int p, int errors, node_t v) const {
        //return false;  // for the intervals to work correctly

        auto crumbs_it = C.find(v); 
        if (crumbs_it != C.end())
            if (crumbs_it->second.contains(p))
                return true;
        return false;
    }

    inline void update_crumbs_for_node(int p, int errors, node_t v) {
        cost_t cost = args.max_seed_errors + 1 - errors;
		auto crumbs_it = C.find(v);
		if (crumbs_it == C.end()) {
			++read_cnt.states_with_crumbs;

			std::unordered_map<int, cost_t> tmp;
			tmp[p] = cost;
			C[v] = tmp;
		}
		else {
			auto &seeds = crumbs_it->second;
			if (seeds.contains(p)) {
				++read_cnt.repeated_states;    
			} else {
				seeds[p] = cost;
				++read_cnt.states_with_crumbs;
			}
		}
		assert(C.contains(v) && C[v].contains(p));
    }

    void update_crumbs_up_the_trie(int p, int errors, node_t trie_v) {
        assert(G.node_in_trie(trie_v));
        ++read_cnt.paths_considered;

        do {
            if (crumbs_already_set(p, errors, trie_v))  // optimization
                return;

            update_crumbs_for_node(p, errors, trie_v);
            if (trie_v == 0)
                break;
            int cnt = 0;
            for (auto it=G.begin_orig_rev_edges(trie_v); it!=G.end_orig_rev_edges(); ++it) {
                trie_v = it->to;
                assert (G.node_in_trie(trie_v));
                ++cnt;
            }
            assert(trie_v == 0 || cnt == 1);
        } while (true);
    }

    // `C[u]++` for all nodes (incl. `v`) that lead from `0` to `v` with a path of length `i`
    // Fully ignores labels.
    // Returns if the the supersource was reached at least once.
    bool put_crumbs_backwards(int p, int i, node_t v, int errors) {
        for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it) {
            if (!crumbs_already_set(p, errors, it->to)) {  // Checking leafs is enough
                if (G.node_in_trie(it->to)) {  // If goint go try -> skip the queue
                    if (i-1 - G.get_trie_depth() <= max_indels_) {
						update_crumbs_up_the_trie(p, errors, it->to);
                    }
                } else if (i-1 >= -max_indels_ + G.get_trie_depth()) { // prev in GRAPH
					update_crumbs_for_node(p, errors, it->to);
                    put_crumbs_backwards(p, i-1, it->to, errors);
                }
            }
        }

        return true;
    }

    // Assumes that seed_len <= D so there are no duplicating outgoing labels.
    void match_seed_put_crumbs(const read_t *r, int p, int start, int i, node_t v, int remaining_errors) {
        if (i < start + args.seed_len) {
            // Match exactly down the trie and then through the original graph.
            for (auto it=G.begin_all_matching_edges(v, r->s[i]); it!=G.end_all_edges(); ++it) {
                int new_remaining_errors = remaining_errors;
                int new_i = i;
                if (it->type != DEL)
                    ++new_i;
                if (it->type != ORIG && it->type != JUMP)  // ORIG in the graph, JUMP in the trie
                    --new_remaining_errors;
                if (new_remaining_errors >= 0)
                    match_seed_put_crumbs(r, p, start, new_i, it->to, new_remaining_errors);
            }
//        } else if (G.node_in_trie(v)) {
//            // Climb the trie.
//            for (auto it=G.begin_orig_edges(v); it!=G.end_orig_edges(); ++it)
//                match_seed_put_crumbs(r, start, i+1, it->to);
        } else {
            assert(!G.node_in_trie(v));
			put_crumbs_backwards(p, i, v, args.max_seed_errors-remaining_errors);
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
            match_seed_put_crumbs(r, seeds, i, i, 0, args.max_seed_errors);  // seed from [i, i+seed_len)
            seeds++;
        }
        return seeds;
    }
};

}
