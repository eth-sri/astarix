#pragma once

#include <map>
#include <string>
#include <vector>

#include "graph.h"
#include "utils.h"
#include "io.h"

typedef int node_t;

namespace astarix {

class AStarLandmarks: public AStarHeuristic {
  private:
    // Fixed parameters
    const graph_t &G;
    const EditCosts &costs;

    // Updated separately for every read
    int pivots, pivot_length;
    const read_t *last_r;

    // H[u] := number of exactly aligned pivots after node `u`
    // It is safe to increase more to H than needed.
    std::vector<int> H;
//    std::vector<int> visited;  // == +1 if a node has already been added to H; -1 if a node has already been subtracted from H
    //std::unordered_map<node_t, cost_t> _star;

    // stats
    int reads;
    int landmark_matches;

  public:
    AStarLandmarks(const graph_t &_G, const EditCosts &_costs)
        : G(_G), costs(_costs) {
        H.resize(G.nodes());
        pivot_length = 30;
        reads = 0;
        landmark_matches = 0;
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
    // O(1)
    cost_t h(const state_t &st) const {
        int cnt = __builtin_popcount(H[st.v]);
        assert(pivots >= cnt);
        LOG_DEBUG << "h(" << st.i << ", " << st.v << ") = ("
                  << pivots << "-" << cnt << ") * " << costs.get_min_mismatch_cost();
        return (pivots - cnt) * costs.get_min_mismatch_cost();
    }

    // Cut r into chunks of length pivot_len, starting from the end.
    void before_every_alignment(const read_t *r) {
        ++reads;
        last_r = r;
        pivots = gen_pivots_and_update(r, pivot_length, +1);
    }

    void after_every_alignment() {
        // Revert the updates by adding -1 instead of +1.
        gen_pivots_and_update(last_r, pivot_length, -1);

        // TODO: removedebug
        for(int i=0; i<H.size(); i++)
            assert(H[i] == 0);
    }

    void print_params(std::ostream &out) const {
    }

    void print_stats(std::ostream &out) const {
        out << " == AStarixLandmarks == " << std::endl;
        out << " Total landmark matches for all reads: " << landmark_matches
            << "(" << 1.0*landmark_matches/reads << " per read)" << std::endl;
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
    // Returns if the the supersource was reached at least once.
    bool update_path_backwards(int p, int i, node_t v, int dval) {
        LOG_DEBUG << "Backwards trace: (" << i << ", " << v << ")";
        if (dval == +1) H[v] |= 1<<p;   // fire p-th bit
        else H[v] &= ~(1<<p);  // remove p-th bit
        LOG_DEBUG << "H[" << v << "] = " << H[v];
        //assert(H[v] >= 0);
        //assert(H[v] <= pivots);

        if (i == 0) {
            assert(v == 0);  // supersource is reached; no need to update H[0]
            return true;
        }

        bool at_least_one_path = false;
        for (auto it=G.begin_orig_rev_edges(v); it!=G.end_orig_rev_edges(); ++it) {
            // TODO: iterate edges from the trie separately from the edges from the graph
            if (should_proceed_backwards_to(i-1, it->to)) {
                assert(update_path_backwards(p, i-1, it->to, dval));
                at_least_one_path = true; // debug
            }
        }

        return at_least_one_path;
    }

    // Assumes that pivot_len <= D so there are no duplicating outgoing labels.
    // In case of success, v is the 
    // Returns a list of (shift_from_start, v) with matches
    void match_pivot_and_update(const read_t *r, int p, int start, int pivot_len, int i, node_t v, int dval) {
        if (i < start + pivot_len) {
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
            LOG_DEBUG << "Updating for pivot (" << i << ", " << v << ") with dval=" << dval;
            assert(update_path_backwards(p, i, v, dval));

            ++landmark_matches;  // debug info
        }
    }

    // Split r into pivots of length pivot_len.
    // For each exact occurence of a pivot (i,v) in the graph,
    //   add dval to H[u] for all nodes u on the path of match-length exactly `i` from supersource `0` to `v`
    // Returns the number of pivots.
    int gen_pivots_and_update(const read_t *r, int pivot_len, int dval) {
        int pivots = 0;
        for (int i=r->len-pivot_len+1; i>0; i-=pivot_len) {
            // pivot from [i, i+pivot_len)
            match_pivot_and_update(r, pivots, i, pivot_len, i, 0, dval);
            pivots++;
        }
        return pivots;
    }
};

}
