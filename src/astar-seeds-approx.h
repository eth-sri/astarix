#pragma once

#include <map>
#include <string>
#include <vector>

#include "graph.h"
#include "utils.h"
#include "io.h"

typedef int node_t;

namespace astarix {

class AStarSeedsWithErrors: public AStarHeuristic {
  private:
    // Fixed parameters
    const graph_t &G;
    const EditCosts &costs;
    const arguments::AStarSeedsArgs args;

    // Updated separately for every read
    const read_t *r_;

    // H[u] := number of exactly aligned seeds after node `u`
    // It is safe to increase more to H than needed.
    // TODO: make H dependent on the distance to the seed
    static const int MAX_SEED_ERRORS = 3;
    std::vector<int> H[MAX_SEED_ERRORS+1];
//    std::vector<int> visited;  // == +1 if a node has already been added to H; -1 if a node has already been subtracted from H
    //std::unordered_map<node_t, cost_t> _star;

    struct Counters {
        Counter seeds;          // number of seeds (depends only on the read)
        Counter seed_matches;   // places in the graph where seeds match well enough
        Counter paths_considered;  // number of paths updated for all seeds (supersource --> match)
        Counter marked_states;     // the number of states marked over all paths for all seeds
        Counter repeated_states;   // number of times crumbs are put on a states that already has crumbs

        void clear() {
            seeds.clear();
            seed_matches.clear();
            paths_considered.clear();
            marked_states.clear();
            repeated_states.clear();
        }

        Counters& operator+=(const Counters &b) {
            seeds += b.seeds;
            seed_matches += b.seed_matches;
            paths_considered += b.paths_considered;
            marked_states += b.marked_states;
            repeated_states += b.repeated_states;
            return *this;
        }
    };

    Counters read_cnt, global_cnt;

    Counter reads_;            // reads processed

    cost_t best_heuristic_sum_;

    // debug; TODO: remove
    std::unordered_set<node_t> visited_nodes_backwards;

  public:
    AStarSeedsWithErrors(const graph_t &_G, const EditCosts &_costs, arguments::AStarSeedsArgs _args)
        : G(_G), costs(_costs), args(_args) {
        for (int i=0; i<=args.max_seed_errors; i++)
            H[i].resize(G.nodes());
        //LOG_INFO << "A* matching class constructed with:";
        //LOG_INFO << "  seed_len    = " << args.seed_len;
    }

    // assume only hamming distance (substitutions)
    // every seed is a a pair (s,j), s.t. s=r[j...j+seed_len)
    //                              j2       j1       j0
    // r divided into seeds: ----|---s2---|---s1---|---s0---|
    // alignment:                u  v2       v1       v0
    //
    // h(<u,i>) := P - f(<u,i>)
    // f(<u,i>) := |{ (s,j) \in seed | exists v: exists path from u->v of length exactly (j-i) and s aligns exactly from v }|,
    // 
    // where P is the number of seeds
    // Accounts only for the last seeds.
    // O(1)
    cost_t h(const state_t &st) const {
        int all_seeds_to_end = (r_->len - st.i - 1) / args.seed_len;

        int total_errors = (args.max_seed_errors+1)*all_seeds_to_end;  // the maximum number of errors
        int not_used_mask = ((1<<all_seeds_to_end)-1);   // at first no seed is used: 111...11111 (in binary)
        for (int errors=0; errors<=args.max_seed_errors; errors++) {
            int h_remaining = H[errors][st.v] & not_used_mask;
            int matching_seeds = __builtin_popcount(h_remaining);
            assert(matching_seeds <= all_seeds_to_end);
            not_used_mask &= ~h_remaining;  // remove the bits for used seeds

            total_errors -= matching_seeds*(args.max_seed_errors+1-errors); // for 0 errors, lower the heuristic by max_seed_errors+1
        }
        
        cost_t res = total_errors * costs.get_min_mismatch_cost();
        //LOG_DEBUG << "h(" << st.i << ", " << st.v << ") = "
        //          << total_errors << " * " << costs.get_min_mismatch_cost()
        //          << " = " << res;
        return res;
    }

    // Cut r into chunks of length seed_len, starting from the end.
    void before_every_alignment(const read_t *r) {
        ++reads_;

        read_cnt.clear();
        read_cnt.seeds.set(gen_seeds_and_update(r, +1));

        r_ = r;

        global_cnt += read_cnt;

        best_heuristic_sum_ += h(state_t(0.0, 0, 0, -1, -1));

        LOG_INFO << r->comment << " A* seeds stats: "
            << read_cnt.seeds << " seeds " 
            << "matching at " << read_cnt.seed_matches << " graph positions "
            << "and generating " << read_cnt.paths_considered << " paths "
            << "over " << read_cnt.marked_states << " states"
            << "(" << read_cnt.repeated_states << " repeated)"
            << "with best heuristic " << h(state_t(0.0, 0, 0, -1, -1)) << " "
            << "out of possible " << (args.max_seed_errors+1)*read_cnt.seeds.get();
//             << "which compensates for <= " << (args.max_seed_errors+1)*read_cnt.seeds - h(state_t(0.0, 0, 0, -1, -1)) << " errors";
    }

    void after_every_alignment() {
        // Revert the updates by adding -1 instead of +1.
        gen_seeds_and_update(r_, -1);

        // TODO: removedebug
        for(int e=0; e<MAX_SEED_ERRORS; e++)
            for(int i=0; i<H[e].size(); i++)
                assert(H[e][i] == 0);
    }

    void print_params(std::ostream &out) const {
        out << "     seed length: " << args.seed_len << " bp"                    << std::endl;
        out << " max seed errors: " << args.max_seed_errors                      << std::endl;
        out << "     shifts allowed: " << args.max_indels                        << std::endl;
    }

    void print_stats(std::ostream &out) const {
        out << "        For all reads:"                                                     << std::endl;
        out << "                            Seeds: " << global_cnt.seeds                    << std::endl;
        out << "                     Seed matches: " << global_cnt.seed_matches
            << "(" << 1.0*global_cnt.seed_matches.get()/reads_.get() << " per read)"                   << std::endl;
        out << "                    Paths considered: " << global_cnt.paths_considered      << std::endl;
        out << "                       States marked: " << global_cnt.marked_states << "(" << global_cnt.repeated_states << " repeated)"  << std::endl;
        out << "                Best heuristic (avg): " << (double)best_heuristic_sum_/reads_.get() << std::endl;
//        out << "                 Max heursitic value: " << ?? << std::endl;
    }

  private:
    // Proceed if node is in the trie and it is time for the trie,
    //   or if the node is not in the trie and it is too early for the trie.
    // trie_depth is the number of node levels (including the supersource)
    //bool should_proceed_backwards_to(int i, node_t v) {
    //    bool time_for_trie = i < G.get_trie_depth();
    //    bool node_in_trie = G.node_in_trie(v);
    //    //LOG_DEBUG << "next node: " << v << ", time_for_trie: " << time_for_trie << ", node_in_trie: " << node_in_trie;
    //    return !(time_for_trie ^ node_in_trie);
    //}

    inline int diff(int a, int b) { return abs(a-b); }

    const static int CYCLE_VAL = -999999;

    inline bool time_for_trie(int min_i, int max_shifts) {
        return min_i == CYCLE_VAL || min_i <= max_shifts + G.get_trie_depth();
    }

    inline bool crumbs_already_set(int p, int dval, int errors, node_t v) {
        return (dval == -1) ^ bool(H[errors][v] & (1<<p));
    }

    inline void update_crumbs_for_node(int p, int dval, int errors, node_t v) {
        if (dval == +1) {
            if (!(H[errors][v] & (1<<p))) {
                ++read_cnt.marked_states;
            } else {
                ++read_cnt.repeated_states;    
            }
            H[errors][v] |= 1<<p;   // fire p-th bit
        } else {
            H[errors][v] &= ~(1<<p);  // remove p-th bit
        }
    }

    void update_crumbs_up_the_trie(int p, int dval, int errors, node_t v) {
        // TODO: stop in case the crumbs is already set
        // from TRIE
        ++read_cnt.paths_considered;
        if (crumbs_already_set(p, dval, errors, v))
            return;

        update_crumbs_for_node(p, dval, errors, v);
        while (v != 0) { // while not in root
            for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it) { // TODO: use up_trie iterator
                if (G.node_in_trie(it->to)) {
                    v = it->to;
                    if (crumbs_already_set(p, dval, errors, v))
                        return;
                    update_crumbs_for_node(p, dval, errors, v);
                    break;
                }
            }
        }
    }

    void propagate_cycle_backwards(std::unordered_map<node_t, int> *min_i, int p, int dval, int errors, int i, node_t v, int max_shifts) {
        auto curr_min_i_it = min_i->find(v);
        if (curr_min_i_it == min_i->end()) {
            (*min_i)[v] = CYCLE_VAL;
            update_crumbs_for_node(p, dval, errors, v);
        }
        else if (curr_min_i_it->second == CYCLE_VAL)
            return;

        // TODO: add trie
        for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it)
            if ( ( G.node_in_trie(it->to) ) // always go up the TRIE
              || (!G.node_in_trie(it->to) && i-1 >= -max_shifts + G.get_trie_depth())) // prev in GRAPH
                propagate_cycle_backwards(min_i, p, dval, errors, i-1, it->to, max_shifts);
    }

    bool update_path_backwards_bfs(int p, int start_i, node_t start_v, int dval, int max_shifts, int errors) {
        LOG_DEBUG_IF(dval == +1) << "Backwards BFS: p=" << p << ", start_i=" << start_i << ", start_v=" << start_v << ", dval=" << dval << ", max_shifts=" << max_shifts << ", errors=" << errors;

        std::queue< std::pair<node_t, int> > Q;  // (v, i)
        std::unordered_map<node_t, int> min_i;  // min_i[u] = min_{(u,v) \in RGE}(min_i[v] - 1) <= i if min_i[u] isn't defined; or CYCLE_VAL otherwise (before a cycle)

        // init BFS
        Q.push( std::make_pair(start_v, start_i) );

        while(!Q.empty()) {
            auto [v, i] = Q.front(); Q.pop();

            if (v == 0) {
                ++read_cnt.paths_considered;
                continue;
            }

            int curr_min_i=-123123123;

            LOG_DEBUG_IF(dval == +1) << "v=" << v << ", i=" << i << ", H[" << errors << "][" << v << "] = " << H[errors][v];

            auto curr_min_i_it = min_i.find(v);
            if (curr_min_i_it == min_i.end())          curr_min_i = i;   // init with i
            else {
                // already existing
                if      (curr_min_i_it->second == CYCLE_VAL)    continue;         // reaching a cycle which has taken care of propagating, TODO: continue!!
                else if (i < curr_min_i_it->second) {
                    propagate_cycle_backwards(&min_i, p, dval, errors, i, v, max_shifts);  // new cycle found -- propage it back
                    curr_min_i = CYCLE_VAL;
                }
                else {
                    assert(i == curr_min_i_it->second);  // BFS does shortest paths first so the only option is to have two paths with the same length
                    continue;         // other path to the same node have higher length because of BFS
                }
            }

            assert(curr_min_i != -123123123);

            LOG_DEBUG_IF(dval == +1) << "curr_min_i=" << curr_min_i;
            min_i[v] = curr_min_i;
    
            //if ( G.node_in_trie(v) ) {
            //    update_crumbs_up_the_trie(p, dval, errors, v);
            //} else {
                update_crumbs_for_node(p, dval, errors, v);
                for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it) {
                    //LOG_DEBUG_IF(dval == +1) << "Traverse the reverse edge " << v << "->" << it->to << " with label " << it->label;
                    if ( ( G.node_in_trie(it->to) && time_for_trie(curr_min_i-1, max_shifts)) // up the TRIE
                      || (!G.node_in_trie(it->to) && i-1 >= -max_shifts + G.get_trie_depth())) { // prev in GRAPH
                        Q.push( std::make_pair(it->to, i-1) );
                    }
                }
            //}
        }

        return true;
    }

    // `H[u]+=dval` for all nodes (incl. `v`) that lead from `0` to `v` with a path of length `i`
    // Fully ignores labels.
    // Returns if the the supersource was reached at least once.
    bool update_path_backwards_dfs(int p, int i, node_t v, int dval, int max_shifts, int errors) {
        LOG_DEBUG_IF(dval == +1) << "Backwards trace: (" << i << ", " << v << ")";
        update_crumbs_for_node(p, dval, errors, v);

        //LOG_DEBUG_IF(dval == +1) << "H[" << errors << "][" << v << "] = " << H[errors][v];
        //assert(__builtin_popcount(H[errors][v]) <= seeds);

        // linear DFS instead of exponentially many paths
        if (visited_nodes_backwards.find(v) != visited_nodes_backwards.end()) {
            return false;
        } else {
            visited_nodes_backwards.insert(v);
        }

        if (v == 0) {
//            assert(i == 0);  // supersource is reached; no need to update H[0]
            ++read_cnt.paths_considered;
            return true;
        }

        bool at_least_one_path = false;
        for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it) {
            // TODO: iterate edges from the trie separately from the edges from the graph
            LOG_DEBUG_IF(dval == +1) << "Traverse the reverse edge " << v << "->" << it->to << " with label " << it->label;
            if ( (G.node_in_trie(v))  // (1) already in trie
              || (diff(i-1, G.get_trie_depth()) <= max_shifts)  // (2) time to go to trie
              || (i-1 > G.get_trie_depth() && !G.node_in_trie(it->to))) {  // (3) proceed back
                bool success = update_path_backwards_dfs(p, i-1, it->to, dval, max_shifts, errors);
                if (success) at_least_one_path = true;
            }
        }

        return at_least_one_path;
    }

    // Assumes that seed_len <= D so there are no duplicating outgoing labels.
    void match_seed_and_update(const read_t *r, int p, int start, int i, node_t v, int dval, int remaining_errors) {
        //LOG_DEBUG_IF(dval == +1) << "Match forward seed " << p << "[" << start << ", " << start+seed_len << ") to state (" << i << ", " << v << ")"
        //                         << " with " << remaining_errors << " remaining errors.";
        if (i < start + args.seed_len) {
            // TODO: match without recursion if unitig
            // Match exactly down the trie and then through the original graph.
            for (auto it=G.begin_all_matching_edges(v, r->s[i]); it!=G.end_all_edges(); ++it) {
                int new_remaining_errors = remaining_errors;
                int new_i = i;
                //LOG_DEBUG << "edge: " << it->label << ", " << edgeType2str(it->type);
                if (it->type != DEL)
                    ++new_i;
                if (it->type != ORIG && it->type != JUMP)  // ORIG in the graph, JUMP in the trie
                    --new_remaining_errors;
                if (new_remaining_errors >= 0)
                    match_seed_and_update(r, p, start, new_i, it->to, dval, new_remaining_errors);
            }
//        } else if (G.node_in_trie(v)) {
//            // Climb the trie.
//            for (auto it=G.begin_orig_edges(v); it!=G.end_orig_edges(); ++it)
//                match_seed_and_update(r, start, i+1, it->to, dval);
        } else {
            assert(!G.node_in_trie(v));
            //LOG_INFO_IF(dval == +1) << "Updating for seed " << p << "(" << i << ", " << v << ") with dval=" << dval << " with " << max_seed_errors-remaining_errors << " errrors.";
            visited_nodes_backwards.clear();
            bool success = update_path_backwards_bfs(p, i, v, dval, args.max_indels, args.max_seed_errors-remaining_errors);
            //bool success = update_path_backwards_dfs(p, i, v, dval, args.max_indels, args.max_seed_errors-remaining_errors);
            assert(success);

            ++read_cnt.seed_matches;  // debug info
        }
    }

    // Split r into seeds of length seed_len.
    // For each exact occurence of a seed (i,v) in the graph,
    //   add dval to H[u] for all nodes u on the path of match-length exactly `i` from supersource `0` to `v`
    // Returns the number of seeds.
    int gen_seeds_and_update(const read_t *r, int dval) {
        int seeds = 0;
        for (int i=r->len-args.seed_len; i>=0; i-=args.seed_len) {
            // seed from [i, i+seed_len)
            match_seed_and_update(r, seeds, i, i, 0, dval, args.max_seed_errors);
            seeds++;
        }
        return seeds;
    }
};

}
