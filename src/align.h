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
#include "utils.h"

#include <plog/Log.h>

namespace astarix {

struct pairhash {
  public:
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const {
        return std::hash<T>()(x.first) ^ (std::hash<U>()(x.second) + 47);
    }
};

struct Stats {
    Counter<> pushed, popped, greedy_matched;
    Counter<> popped_trie, popped_ref;
    Counter<> explored_states;
    Counter<> repeated_visits;
    AlignerTimers t;

    struct AlignStatus {
        Counter<cost_t> cost;
        Counter<> unique, ambiguous, overcost;

        void clear() {
            cost.clear();
            unique.clear();
            ambiguous.clear();
            overcost.clear();
        }
        int total() const {
            return unique.get() + ambiguous.get() + overcost.get();
        }
        int aligned() const {
            return unique.get() + ambiguous.get();
        }
        AlignStatus& operator+=(const AlignStatus &b) {
            if (b.aligned())
                cost += b.cost;
            unique += b.unique;
            ambiguous += b.ambiguous;
            overcost += b.overcost;
            return *this;
        }
    } align_status;

    //std::vector<Counter> pushed_hist;  // pushed_hist[i] -- number of pushed states <*,i>

    void clear() {
        pushed.clear();
        popped.clear();
        greedy_matched.clear();
        popped_trie.clear();
        popped_ref.clear();
        explored_states.clear();
        repeated_visits.clear();
        align_status.clear();
        t.clear();

    //    pushed_hist.clear();
    }

    Stats& operator+=(const Stats &b) {
        pushed += b.pushed;
        popped += b.popped;
        greedy_matched += b.greedy_matched;
        popped_trie += b.popped_trie;
        popped_ref += b.popped_ref;
        explored_states += b.explored_states;
        repeated_visits += b.repeated_visits;
        align_status += b.align_status;
        t += b.t;

     //   if (pushed_hist.size() < b.pushed_hist.size())
     //       pushed_hist.resize(b.pushed_hist.size());
     //   for (int i=0; i<pushed_hist.size(); i++)
     //       if (i < b.pushed_hist.size())
     //           pushed_hist[i] += b.pushed_hist[i];

        return *this;
    }
};

struct AlignParams {
    const EditCosts &costs;
    const bool greedy_match;
    const cost_t max_align_cost;

    AlignParams(const EditCosts &_costs, const bool _fast_forward, const cost_t _max_align_cost)
      : costs(_costs),
        greedy_match(_fast_forward),
        max_align_cost(_max_align_cost) {
    }

    void print() const {
        LOG_INFO << "Params: ";
        LOG_INFO << "  greedy_match  = " << greedy_match;
        LOG_INFO << "Edit costs: ";
        LOG_INFO << "  match_cost    = " << (int)costs.match;
        LOG_INFO << "  mismatch_cost = " << (int)costs.subst;
        LOG_INFO << "  ins_cost      = " << (int)costs.ins;
        LOG_INFO << "  del_cost      = " << (int)costs.del;
        LOG_INFO << "Max align cost  = " << (int)max_align_cost;
    }
};

class Aligner {
    const graph_t &G;
    const AlignParams &params;

    typedef std::unordered_map<std::pair<pos_t,node_t>, state_t, pairhash> path_t;
    typedef std::unordered_map<std::pair<pos_t,node_t>, edge_t, pairhash> prev_edge_t;
    typedef std::unordered_map<std::pair<pos_t,node_t>, bool, pairhash> visited_t;

  public:
    // Local vars
    path_t p;
    prev_edge_t pe;
	visited_t vis;

    AStarHeuristic *astar;    // Concurrent Aligner's can read and write to the same AStar (it computes and memoizes heuristics).

    mutable Stats stats;

    Aligner(const graph_t &_G, const AlignParams &_params, AStarHeuristic *_astar)
            : G(_G), params(_params), astar(_astar) {
    }

    inline const graph_t& graph() const {
        return G;
    }

    inline const AStarHeuristic& get_astar() const {
        return *astar;
    }

    inline void astar_before_every_alignment(const read_t *r) {
		// Clean up
		p.clear();
		pe.clear();
		vis.clear();

		stats.clear();
		stats.t.total.start();
		//stats.pushed_hist.resize(r.len);

        stats.t.astar_prepare_reads.start();
        astar->before_every_alignment(r);
        stats.t.astar_prepare_reads.stop();
    }

    inline void astar_after_every_alignment() {
        assert(stats.align_status.total() == 1);

		auto timers_copy = stats.t;  // in order not to pass a running timer as an argument
		timers_copy.total.stop();

        stats.t.astar_prepare_reads.start();
        astar->after_every_alignment(timers_copy);
        stats.t.astar_prepare_reads.stop();

		stats.t.total.stop();
    }
  
  private:
    inline void push(queue_t &Q, cost_t sort_cost, const state_t &st) {
        stats.pushed.inc();
        stats.explored_states.inc();
        //stats.pushed_hist[st.i].inc();
		Q.push(score_state_t(sort_cost, st));
    }

    inline score_state_t pop(queue_t &Q) {
        stats.popped.inc();
        score_state_t el = Q.top();
		Q.pop();
        return el;
    }

    inline const state_t& get_const_path(const path_t &p, pos_t i, node_t v) const {
        auto it = p.find(std::make_pair(i, v));
        assert(it != p.end());
        return it->second;
    }

    inline state_t& get_path(path_t &p, pos_t i, node_t v) {
        auto &ret = p[std::make_pair(i, v)];
        return ret;
    }

    inline const edge_t& get_prev_edge(const prev_edge_t &pe, pos_t i, node_t v) const {
        auto it = pe.find(std::make_pair(i, v));
        assert(it != pe.end());
        return it->second;
    }

    inline void set_prev_edge(prev_edge_t &pe, pos_t i, node_t v, const edge_t &e) {
        pe[std::make_pair(i, v)] = e;
    }

    inline bool& visited(visited_t &vis, pos_t i, node_t v) {
        auto &ret = vis[std::make_pair(i, v)];
        return ret;
    }

  public:
    void get_best_path_to_state(state_t final_state, edge_path_t *best_path) const {
        assert(best_path->empty());

        for (state_t curr=final_state; curr.prev_i != -1 && curr.prev_v != -1; curr = get_const_path(p, curr.prev_i, curr.prev_v)) {
            edge_t e = get_prev_edge(pe, curr.i, curr.v);
            best_path->push_back(e);
        }
        reverse(best_path->begin(), best_path->end());
    }

    state_t proceed_identity(path_t &p, prev_edge_t &pe, state_t curr, const read_t &r);

    void try_edge(const read_t &r, const state_t &curr, path_t &p, prev_edge_t &pe, const std::string &algo, queue_t &Q, const edge_t &e);

    /*** A-star and Dijkstra logic ***
        f(n) = g(n) + h(n)
        f(n) -- defines the order in the queue, i.e. the order to pop
        g(n) -- the optimal length from the start to `n'
        h(n) -- a lower bound on the length from `n' to the end
    
        r is a 1-based query
    */
    std::vector<state_t> readmap(const read_t &r, std::string algo, int max_best_alignments);
};

}
