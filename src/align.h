#pragma once

#include <algorithm>
#include <vector>
#include <cstring>
#include <string>
#include <queue>
#include <map>
//#include <sparsehash/sparse_hash_map>
#include <unordered_map>

#include "graph.h"

#include <plog/Log.h>

namespace astarix {

struct pairhash {
  public:
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const {
        return std::hash<T>()(x.first) ^ (std::hash<U>()(x.second) + 47);
    }
};

struct Counters {
    Counter pushed, popped, greedy_matched;
    Counter popped_trie, popped_ref;
    Counter explored_states;
    Counter repeated_visits;

    void clear() {
        pushed.clear();
        popped.clear();
        greedy_matched.clear();
        popped_trie.clear();
        popped_ref.clear();
        explored_states.clear();
        repeated_visits.clear();
    }

    Counters& operator+=(const Counters &b) {
        pushed += b.pushed;
        popped += b.popped;
        greedy_matched += b.greedy_matched;
        popped_trie += b.popped_trie;
        popped_ref += b.popped_ref;
        explored_states += b.explored_states;
        repeated_visits += b.repeated_visits;
        return *this;
    }
};

struct AlignerTimers {
    Timer queue, ff, dicts, astar, total;
    Timer astar_prepare_reads;

    void clear() {
        queue.clear();
        ff.clear();
        dicts.clear();
        astar.clear();
        astar_prepare_reads.clear();
        total.clear();
    }

    AlignerTimers& operator+=(const AlignerTimers &b) {
        queue += b.queue;
        ff += b.ff;
        dicts += b.dicts;
        astar += b.astar;
        astar_prepare_reads += b.astar_prepare_reads;
        total += b.total;
        return *this;
    }
};

struct AlignParams {
    const EditCosts &costs;
    const bool greedy_match;

    AlignParams(const EditCosts &_costs, bool _fast_forward)
      : costs(_costs),
        greedy_match(_fast_forward) {
    }

    void print() const {
        LOG_INFO << "Params: ";
        LOG_INFO << "  greedy_match  = " << greedy_match;
        //LOG_INFO << "  tree depth    = " << _tree_depth;
        LOG_INFO << "Edit costs: ";
        LOG_INFO << "  match_cost    = " << (int)costs.match;
        LOG_INFO << "  mismatch_cost = " << (int)costs.subst;
        LOG_INFO << "  ins_cost      = " << (int)costs.ins;
        LOG_INFO << "  del_cost      = " << (int)costs.del;
    }
};

class Aligner {
    const graph_t &G;
    const AlignParams &params;

    typedef std::unordered_map<std::pair<int,int>, state_t, pairhash> path_t;
    typedef std::unordered_map<std::pair<int,int>, edge_t, pairhash> prev_edge_t;
    typedef std::unordered_map<std::pair<int,int>, bool, pairhash> visited_t;       // not needed; for performance analysis only

  public:
    int unique_best;
    AStarHeuristic *astar;    // Concurrent Aligner's can read and write to the same AStar (it computes and memoizes heuristics).

    // total_counters are aggregating read_counters at the end of every alignment
    mutable Counters read_counters;
    mutable AlignerTimers read_timers;

    Aligner(const graph_t &_G, const AlignParams &_params, AStarHeuristic *_astar)
            : G(_G), params(_params), astar(_astar), unique_best(-1) {
    }

    inline const graph_t& graph() const {
        return G;
    }

    inline const AStarHeuristic& get_astar() const {
        return *astar;
    }

    inline void astar_before_every_alignment(const read_t *r) {
        read_timers.astar_prepare_reads.start();
        astar->before_every_alignment(r);
        read_timers.astar_prepare_reads.stop();
    }

    inline void astar_after_every_alignment() {
        read_timers.astar_prepare_reads.start();
        astar->after_every_alignment();
        read_timers.astar_prepare_reads.stop();
    }
  
  private:
    inline void push(queue_t &Q, cost_t sort_cost, const state_t &st) {
        read_counters.pushed.inc();
        read_timers.queue.start();
            Q.push(score_state_t(sort_cost, st));
        read_timers.queue.stop();
    }

    inline state_t pop(queue_t &Q) {
        read_counters.popped.inc();
        read_timers.queue.start();
        score_state_t el = Q.top();
        Q.pop();
        read_timers.queue.stop();
        return el.second;
    }

    inline const state_t& get_const_path(const path_t &p, int i, int v) const {
        read_timers.dicts.start();
        auto it = p.find(std::make_pair(i, v));
        read_timers.dicts.stop();
        assert(it != p.end());
        return it->second;
    }

    inline state_t& get_path(path_t &p, int i, int v) {
        read_timers.dicts.start();
        auto &ret = p[std::make_pair(i, v)];
        read_timers.dicts.stop();
        return ret;
    }

    inline const edge_t& get_prev_edge(const prev_edge_t &pe, int i, int v) const {
        read_timers.dicts.start();
        auto it = pe.find(std::make_pair(i, v));
        read_timers.dicts.stop();
        assert(it != pe.end());
        return it->second;
    }

    inline void set_prev_edge(prev_edge_t &pe, int i, int v, const edge_t &e) {
        read_timers.dicts.start();
        pe[std::make_pair(i, v)] = e;
        read_timers.dicts.stop();
    }

    inline bool& visited(visited_t &vis, int i, int v) {
        read_timers.dicts.start();
        auto &ret = vis[std::make_pair(i, v)];
        read_timers.dicts.stop();
        return ret;
    }

    void get_best_path_to_state(const path_t &p, const prev_edge_t &pe, state_t final_state, edge_path_t *best_path) const {
        assert(best_path->empty());

        for (state_t curr=final_state; curr.prev_i != -1 && curr.prev_v != -1; curr = get_const_path(p, curr.prev_i, curr.prev_v)) {
            edge_t e = get_prev_edge(pe, curr.i, curr.v);
            best_path->push_back(e);
        }
        reverse(best_path->begin(), best_path->end());
    }

  public:
    state_t proceed_identity(path_t &p, prev_edge_t &pe, state_t curr, const read_t &r);

    void try_edge(const read_t &r, const state_t &curr, path_t &p, prev_edge_t &pe, const std::string &algo, queue_t &Q, const edge_t &e);

    /*** A-star and Dijkstra logic ***
        f(n) = g(n) + h(n)
        f(n) -- defines the order in the queue, i.e. the order to pop
        g(n) -- the optimal length from the start to `n'
        h(n) -- a lower bound on the length from `n' to the end
    
        r is a 1-based query
    */
    state_t readmap(const read_t &r, std::string algo, edge_path_t *best_path);
};

}
