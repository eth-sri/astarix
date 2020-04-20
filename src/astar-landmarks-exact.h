#pragma once

#include <map>
#include <string>
#include <vector>

#include "graph.h"
#include "utils.h"
#include "io.h"

typedef int node_t;

namespace astarix {

///////////////////////////////////////
// in case max_waymark_errors == 0
///////////////////////////////////////

class AStarLandmarksExact: public AStarHeuristic {
  private:
    // Fixed parameters
    const graph_t &G;
    const EditCosts &costs;

    // Updated separately for every read
    int pivot_len;
    const read_t *r_;
    int shifts_allowed_;  // number of deletions accomodated in first trie_depth nucleotides

    // H[u] := number of exactly aligned pivots after node `u`
    // It is safe to increase more to H than needed.
    // TODO: make H dependent on the distance to the landmark
    std::vector<int> H;
//    std::vector<int> visited;  // == +1 if a node has already been added to H; -1 if a node has already been subtracted from H
    //std::unordered_map<node_t, cost_t> _star;

    // stats
    int waymarks;
    int waymark_matches;
    int paths_considered;
    int marked_states;

    // stats for all reads
    int reads;
    int waymarks_;
    int waymark_matches_;
    int paths_considered_;
    int marked_states_;

  public:
    AStarLandmarksExact(const graph_t &_G, const EditCosts &_costs, int _pivot_len, int _shifts_allowed)
        : G(_G), costs(_costs), pivot_len(_pivot_len), shifts_allowed_(_shifts_allowed)  {
        H.resize(G.nodes());
        reads = 0;
        waymarks_ = 0;
        waymark_matches_ = 0;
        paths_considered_ = 0;
        marked_states_ = 0;
        //LOG_INFO << "A* matching class constructed with:";
        //LOG_INFO << "  pivot_len    = " << pivot_len;
    }

    // assume only hamming distance (substitutions)
    // every landmark is a a pair (s,j), s.t. s=r[j...j+pivot_len)
    //                            j2       j1       j0
    // r divided into pivots: ----|---s2---|---s1---|---s0---|
    // alignment:              u  v2       v1       v0
    //
    // h(<u,i>) := P - f(<u,i>)
    // f(<u,i>) := |{ (s,j) \in pivot | exists v: exists path from u->v of length exactly (j-i) and s aligns exactly from v }|,
    // 
    // where P is the number of pivots
    // Accounts only for the last landmarks.
    // O(1)
    cost_t h(const state_t &st) const {
        int max_landmarks_to_end = (r_->len - st.i - 1) / pivot_len;
        int h_remaining = H[st.v] & ((1<<max_landmarks_to_end)-1);
        int matching_landmarks = __builtin_popcount(h_remaining);
        assert(max_landmarks_to_end >= matching_landmarks);
        cost_t res = (max_landmarks_to_end - matching_landmarks) * costs.get_min_mismatch_cost();
        LOG_DEBUG << "h(" << st.i << ", " << st.v << ") = ("
                  << max_landmarks_to_end << "-" << matching_landmarks << ") * "
                  << costs.get_min_mismatch_cost()
                  << " = " << res;
        return res;
    }

    // Cut r into chunks of length pivot_len, starting from the end.
    void before_every_alignment(const read_t *r) {
        ++reads;
        waymark_matches = 0;
        paths_considered = 0;
        marked_states = 0;
        r_ = r;
        waymarks = gen_pivots_and_update(r, pivot_len, +1);

        waymarks_ += waymarks;
        waymark_matches_ += waymark_matches;
        paths_considered_ += paths_considered;
        marked_states_ += marked_states;
    }

    void after_every_alignment() {
        // Revert the updates by adding -1 instead of +1.
        gen_pivots_and_update(r_, pivot_len, -1);

        // TODO: removedebug
        for(int i=0; i<H.size(); i++)
            assert(H[i] == 0);
    }

    void print_params(std::ostream &out) const {
        out << "   landmark length: " << pivot_len << " bp" << std::endl;
    }

    void print_stats(std::ostream &out) const {
        out << "        For all reads:"                                             << std::endl;
        out << "                             Waymarks: " << waymarks_               << std::endl;
        out << " Total landmark matches for all reads: " << waymark_matches_
            << "(" << 1.0*waymark_matches_/reads << " per read)"                    << std::endl;
        out << "                     Paths considered: " << paths_considered_       << std::endl;
        out << "                   Graph nodes marked: " << marked_states_          << std::endl;
    }

  private:
    // Proceed if node is in the trie and it is time for the trie,
    //   or if the node is not in the trie and it is too early for the trie.
    // trie_depth is the number of node levels (including the supersource)
    bool should_proceed_backwards_to(int i, node_t v) {
        bool time_for_trie = i < G.get_trie_depth();
        bool node_in_trie = G.node_in_trie(v);
        LOG_DEBUG << "next node: " << v << ", time_for_trie: " << time_for_trie << ", node_in_trie: " << node_in_trie;
        return !(time_for_trie ^ node_in_trie);
    }

    // `H[u]+=dval` for all nodes (incl. `v`) that lead from `0` to `v` with a path of length `i`
    // TODO: optimize with string nodes
    // Fully ignores labels.
    // Returns if the the supersource was reached at least once.
    bool update_path_backwards(int p, int i, node_t v, int dval, int shifts_remaining) {
        LOG_DEBUG_IF(dval == +1) << "Backwards trace: (" << i << ", " << v << ")";
        if (dval == +1) {
            if (!(H[v] & (1<<p))) {
                ++marked_states;
                H[v] |= 1<<p;   // fire p-th bit
            }
        } else {
            H[v] &= ~(1<<p);  // remove p-th bit
        }
        LOG_DEBUG_IF(dval == +1) << "H[" << v << "] = " << H[v];
        //assert(__builtin_popcount(H[v]) <= waymarks);

        if (i == 0) {
            assert(v == 0);  // supersource is reached; no need to update H[0]
            ++paths_considered;
            return true;
        }

        bool at_least_one_path = false;
        for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it) {
            // TODO: iterate edges from the trie separately from the edges from the graph
            if (should_proceed_backwards_to(i-1, it->to)) {
                // Go back freely to accomodate deletions.
                if (i-1 == G.get_trie_depth() && !G.node_in_trie(it->to)) {
                    if (shifts_remaining > 0)
                        update_path_backwards(p, i, it->to, dval, shifts_remaining-1);
                }

                LOG_DEBUG_IF(dval == +1) << "Traverse the reverse edge " << v << "->" << it->to << " with label " << it->label;
                bool success = update_path_backwards(p, i-1, it->to, dval, shifts_remaining);
                assert(success);
                at_least_one_path = true; // debug
            }
        }

        return at_least_one_path;
    }

    // Assumes that pivot_len <= D so there are no duplicating outgoing labels.
    // In case of success, v is the 
    // Returns a list of (shift_from_start, v) with matches
    void match_pivot_and_update(const read_t *r, int p, int start, int pivot_len, int i, node_t v, int dval) {
        LOG_DEBUG_IF(dval == +1) << "Match forward pivot " << p << "[" << start << ", " << start+pivot_len << ") to state (" << i << ", " << v << ")";
        if (i < start + pivot_len) {
            // TODO: match without recursion if unitig
            // Match exactly down the trie and then through the original graph.
            for (auto it=G.begin_orig_edges(v); it!=G.end_orig_edges(); ++it)
                if (it->label == r->s[i])
                    match_pivot_and_update(r, p, start, pivot_len, i+1, it->to, dval);
//        } else if (G.node_in_trie(v)) {
//            // Climb the trie.
//            for (auto it=G.begin_orig_edges(v); it!=G.end_orig_edges(); ++it)
//                match_pivot_and_update(r, start, pivot_len, i+1, it->to, dval);
        } else {
            assert(!G.node_in_trie(v));
            LOG_DEBUG_IF(dval == +1) << "Updating for pivot " << p << "(" << i << ", " << v << ") with dval=" << dval;
            bool success = update_path_backwards(p, i, v, dval, shifts_allowed_);
            assert(success);

            ++waymark_matches;
        }
    }

    // Split r into pivots of length pivot_len.
    // For each exact occurence of a pivot (i,v) in the graph,
    //   add dval to H[u] for all nodes u on the path of match-length exactly `i` from supersource `0` to `v`
    // Returns the number of pivots.
    int gen_pivots_and_update(const read_t *r, int pivot_len, int dval) {
        int waymarks = 0;
        for (int i=r->len-pivot_len; i>=0; i-=pivot_len) {
            // pivot from [i, i+pivot_len)
            match_pivot_and_update(r, waymarks, i, pivot_len, i, 0, dval);
            waymarks++;
        }
        return waymarks;
    }
};

}
