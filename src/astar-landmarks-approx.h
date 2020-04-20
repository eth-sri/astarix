#pragma once

#include <map>
#include <string>
#include <vector>

#include "graph.h"
#include "utils.h"
#include "io.h"

typedef int node_t;

namespace astarix {

class AStarWaymarksWithErrors: public AStarHeuristic {
  private:
    // Fixed parameters
    const graph_t &G;
    const EditCosts &costs;
    const int max_waymark_errors;

    // Updated separately for every read
    int pivot_len;
    const read_t *r_;
    const int shifts_allowed_;  // number of deletions accomodated in first trie_depth nucleotides

    // H[u] := number of exactly aligned waymarks after node `u`
    // It is safe to increase more to H than needed.
    // TODO: make H dependent on the distance to the waymark
    static const int MAX_WAYMARK_ERRORS = 3;
    std::vector<int> H[MAX_WAYMARK_ERRORS+1];
//    std::vector<int> visited;  // == +1 if a node has already been added to H; -1 if a node has already been subtracted from H
    //std::unordered_map<node_t, cost_t> _star;

    // stats per read
    int waymarks;          // number of waymarks (depends only on the read)
    int waymark_matches;   // places in the graph where waymarks match well enough
    int paths_considered;  // number of paths updated for all waymarks (supersource --> match)
    int marked_states;     // the number of states marked over all paths for all waymarks

    // global stats
    int reads_;            // reads processed
    int waymarks_;          
    int waymark_matches_;   
    int paths_considered_;  
    int marked_states_;     
    cost_t best_heuristic_sum_;

  public:
    AStarWaymarksWithErrors(const graph_t &_G, const EditCosts &_costs, int _pivot_len, int _max_waymark_errors, int _shifts_allowed)
        : G(_G), costs(_costs), pivot_len(_pivot_len), max_waymark_errors(_max_waymark_errors), shifts_allowed_(_shifts_allowed) {
        for (int i=0; i<=max_waymark_errors; i++)
            H[i].resize(G.nodes());
        reads_ = 0;
        waymarks_ = 0;
        waymark_matches_ = 0;
        paths_considered_ = 0;
        marked_states_ = 0;
        best_heuristic_sum_ = 0;
        //LOG_INFO << "A* matching class constructed with:";
        //LOG_INFO << "  pivot_len    = " << pivot_len;
    }

    // assume only hamming distance (substitutions)
    // every waymark is a a pair (s,j), s.t. s=r[j...j+pivot_len)
    //                              j2       j1       j0
    // r divided into waymarks: ----|---s2---|---s1---|---s0---|
    // alignment:                u  v2       v1       v0
    //
    // h(<u,i>) := P - f(<u,i>)
    // f(<u,i>) := |{ (s,j) \in pivot | exists v: exists path from u->v of length exactly (j-i) and s aligns exactly from v }|,
    // 
    // where P is the number of waymarks
    // Accounts only for the last waymarks.
    // O(1)
    cost_t h(const state_t &st) const {
        int all_waymarks_to_end = (r_->len - st.i - 1) / pivot_len;

        int total_errors = (max_waymark_errors+1)*all_waymarks_to_end;  // the maximum number of errors
        int not_used_mask = ((1<<all_waymarks_to_end)-1);   // at first no waymark is used: 111...11111 (in binary)
        for (int errors=0; errors<=max_waymark_errors; errors++) {
            int h_remaining = H[errors][st.v] & not_used_mask;
            int matching_waymarks = __builtin_popcount(h_remaining);
            assert(matching_waymarks <= all_waymarks_to_end);
            not_used_mask &= ~h_remaining;  // remove the bits for used waymarks

            total_errors -= matching_waymarks*(max_waymark_errors+1-errors); // for 0 errors, lower the heuristic by max_waymark_errors+1
        }
        
        cost_t res = total_errors * costs.get_min_mismatch_cost();
        //LOG_DEBUG << "h(" << st.i << ", " << st.v << ") = "
        //          << total_errors << " * " << costs.get_min_mismatch_cost()
        //          << " = " << res;
        return res;
    }

    // Cut r into chunks of length pivot_len, starting from the end.
    void before_every_alignment(const read_t *r) {
        ++reads_;
        paths_considered = 0;
        marked_states = 0;
        waymark_matches = 0;

        r_ = r;
        waymarks = gen_waymarks_and_update(r, pivot_len, +1);

        waymarks_ += waymarks;
        waymark_matches_ += waymark_matches;
        paths_considered_ += paths_considered;
        marked_states_ += marked_states;
        best_heuristic_sum_ += h(state_t(0.0, 0, 0, -1, -1));

        LOG_INFO << r->comment << " A* waymarks stats: "
                 << waymarks << " waymarks " 
                 << "matching at " << waymark_matches << " graph positions "
                 << "and generating " << paths_considered << " paths "
                 << "over " << marked_states << " states"
                 << "with best heuristic " << h(state_t(0.0, 0, 0, -1, -1)) << " "
                 << "out of possible " << (max_waymark_errors+1)*waymarks;
//                 << "which compensates for <= " << (max_waymark_errors+1)*waymarks - h(state_t(0.0, 0, 0, -1, -1)) << " errors";
    }

    void after_every_alignment() {
        // Revert the updates by adding -1 instead of +1.
        gen_waymarks_and_update(r_, pivot_len, -1);

        // TODO: removedebug
        for(int e=0; e<MAX_WAYMARK_ERRORS; e++)
            for(int i=0; i<H[e].size(); i++)
                assert(H[e][i] == 0);
    }

    void print_params(std::ostream &out) const {
        out << "     waymark length: " << pivot_len << " bp"                        << std::endl;
        out << " max waymark errors: " << max_waymark_errors                        << std::endl;
        out << "     shifts allowed: " << shifts_allowed_                           << std::endl;
    }

    void print_stats(std::ostream &out) const {
        out << "        For all reads:"                                             << std::endl;
        out << "                            Waymarks: " << waymarks_                << std::endl;
        out << "                     Waymark matches: " << waymark_matches_
            << "(" << 1.0*waymark_matches_/reads_ << " per read)"                   << std::endl;
        out << "                    Paths considered: " << paths_considered_        << std::endl;
        out << "                  Graph nodes marked: " << marked_states_           << std::endl;
        out << "                Best heuristic (avg): " << (double)best_heuristic_sum_/reads_ << std::endl;
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

    int diff(int a, int b) { return abs(a-b); }

    // `H[u]+=dval` for all nodes (incl. `v`) that lead from `0` to `v` with a path of length `i`
    // TODO: optimize with string nodes
    // Fully ignores labels.
    // Returns if the the supersource was reached at least once.
    bool update_path_backwards(int p, int i, node_t v, int dval, int shifts_remaining, int errors) {
        LOG_DEBUG_IF(dval == +1) << "Backwards trace: (" << i << ", " << v << ")";
        if (dval == +1) {
            if (!(H[errors][v] & (1<<p))) {
                ++marked_states;
                H[errors][v] |= 1<<p;   // fire p-th bit
            }
        } else {
            H[errors][v] &= ~(1<<p);  // remove p-th bit
        }
        //LOG_DEBUG_IF(dval == +1) << "H[" << errors << "][" << v << "] = " << H[errors][v];
        //assert(__builtin_popcount(H[errors][v]) <= waymarks);

        if (v == 0) {
//            assert(i == 0);  // supersource is reached; no need to update H[0]
            ++paths_considered;
            return true;
        }

        bool at_least_one_path = false;
        for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it) {
            // TODO: iterate edges from the trie separately from the edges from the graph
            LOG_DEBUG_IF(dval == +1) << "Traverse the reverse edge " << v << "->" << it->to << " with label " << it->label;
            if ( (G.node_in_trie(v))  // (1) already in trie
              || (diff(i-1, G.get_trie_depth()) <= shifts_remaining)  // (2) time to go to trie
              || (i-1 > G.get_trie_depth() && !G.node_in_trie(it->to))) {  // (3) proceed back
                bool success = update_path_backwards(p, i-1, it->to, dval, shifts_remaining, errors);
                if (success) at_least_one_path = true;
            }
        }

        return at_least_one_path;
    }

    // Assumes that pivot_len <= D so there are no duplicating outgoing labels.
    // In case of success, v is the 
    // Returns a list of (shift_from_start, v) with matches
    void match_waymark_and_update(const read_t *r, int p, int start, int pivot_len, int i, node_t v, int dval, int remaining_errors) {
        //LOG_DEBUG_IF(dval == +1) << "Match forward pivot " << p << "[" << start << ", " << start+pivot_len << ") to state (" << i << ", " << v << ")"
        //                         << " with " << remaining_errors << " remaining errors.";
        if (i < start + pivot_len) {
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
                    match_waymark_and_update(r, p, start, pivot_len, new_i, it->to, dval, new_remaining_errors);
            }
//        } else if (G.node_in_trie(v)) {
//            // Climb the trie.
//            for (auto it=G.begin_orig_edges(v); it!=G.end_orig_edges(); ++it)
//                match_waymark_and_update(r, start, pivot_len, i+1, it->to, dval);
        } else {
            assert(!G.node_in_trie(v));
            //LOG_INFO_IF(dval == +1) << "Updating for waymark " << p << "(" << i << ", " << v << ") with dval=" << dval << " with " << max_waymark_errors-remaining_errors << " errrors.";
            bool success = update_path_backwards(p, i, v, dval, shifts_allowed_, max_waymark_errors-remaining_errors);
            assert(success);

            ++waymark_matches;  // debug info
        }
    }

    // Split r into pivots of length pivot_len.
    // For each exact occurence of a pivot (i,v) in the graph,
    //   add dval to H[u] for all nodes u on the path of match-length exactly `i` from supersource `0` to `v`
    // Returns the number of pivots.
    int gen_waymarks_and_update(const read_t *r, int pivot_len, int dval) {
        int waymarks = 0;
        for (int i=r->len-pivot_len; i>=0; i-=pivot_len) {
            // pivot from [i, i+pivot_len)
            match_waymark_and_update(r, waymarks, i, pivot_len, i, 0, dval, max_waymark_errors);
            waymarks++;
        }
        return waymarks;
    }
};

}
