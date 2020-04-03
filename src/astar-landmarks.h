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
    // Parameters
    const graph_t &G;
    const EditCosts &costs;

    // Main data structure
    unordered_map<node_t, cost_t> _star;

    // Stats
    Timer astar_time;

  public:
    AStarLandmarks(const graph_t &_G, const EditCosts &_costs)
        : G(_G), costs(_costs) {
        //LOG_INFO << "A* matching class constructed with:";
        //LOG_INFO << "  pivot_len    = " << pivot_len;

        //void prepare(const read_t &r) {
        int pivot_len = 30;
        int pivots = 0;
        _star.clear();
        for (int i=r.len-pivot_len+1; i>0; i-=pivot_len)
            add_pivot(r, i, pivot_len);
    }

    cost_t h(const read_t &r, const state_t &st) {
    }

    void print_params(std::ostream &out) {
    }

    void print_stats(std::ostream &out) {
    }

  private:
    // Return next vertex.
    node_t match_label(node_t v, label_t l) {
        for (auto it=G.begin_orig_edges(v); it!=G.end_orig_edges(); ++it)
            if (it->label == l)
                return it->to;
        return -1;
    }

    // Assumes that pivot_len <= D so there are no duplicating outgoing labels.
    // In case of success, v is the 
    void exact_match(const read_t &r, int start, int pivot_len, int i, node_t v) {
        if (i == start + pivot_len) {
            go_backwards(v);
        }

        node_t v = 0;  // supersource
        for (int i=start; i<start+pivot_len; i++)
            if ((v = match_label(v, r.s[i])) == -1)
                return false;
        return true;
    }

    void add_pivot(const read_t &r, int start, int pivot_len) {
        exact_match(r, start, pivot_len, start, 0);
    }

    size_t get_pivots() {
        return pivots;
    }
};

}
