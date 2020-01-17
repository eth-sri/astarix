#ifndef ASTARIX_GRAPH_H
#define ASTARIX_GRAPH_H

#include <climits>
#include <fstream>
#include <memory.h>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <plog/Log.h>

#include "utils.h"

namespace astarix {

const cost_t INF  = 1000000000.0;

class state_t {
  public:
    cost_t cost;  // cost of a path
    int i;    // index in the query string
    int v;    // vertex in the graph
	int prev_i, prev_v;

    state_t() : cost(INF), i(-1), v(-1), prev_i(-1), prev_v(-1) {}
    state_t(cost_t _cost, int _i, int _v, int _prev_i, int _prev_v) 
        : cost(_cost), i(_i), v(_v), prev_i(_prev_i), prev_v(_prev_v) {}

	bool undef() {
		return i == -1;
	}

    bool operator<(const state_t &other) const {
        return cost > other.cost;
    }

    bool optimize(const state_t &cand) {
		assert(cand.cost != INF);
        if (cand.cost < cost) {
            cost = cand.cost;
            i = cand.i;
            v = cand.v;
			prev_i = cand.prev_i;
			prev_v = cand.prev_v;
            return true;
        }
        return false;
    }
};

static std::ostream& operator<<(std::ostream& os, const state_t &st) {
	os << "(i=" << st.i << ", v=" << st.v << ", cost=" << st.cost << ")";
	return os;
}

typedef std::vector<state_t> path_t;
typedef std::vector<edge_t> edge_path_t;

struct graph_t {
  //private:
	//mutable cost_t _min_edge_cost;
	//mutable cost_t _min_edit_cost;

  public:
    std::vector<edge_t> E;  // n linked lists emulated in a stack E
    std::vector<int> V;     // E[ V[i] ] -- last added outgoing edge from vertex i
	// if a node with number 0 exists, it is a supersource

	int orig_nodes, orig_edges;

	const char *EdgeTypeStr[5];
	int trie_first_node, trie_nodes, trie_edges;

	graph_t(bool _with_reverse_edges=0)
			: trie_nodes(0), trie_edges(0), orig_nodes(0), orig_edges(0)
			//: with_reverse_edges(_with_reverse_edges)
			{
		V.resize(1, -1);  // 0 preserved for a supersource
		//_min_edge_cost = -1;
		//_min_edit_cost = -1;

		EdgeTypeStr[ORIG] = "ORIG";
		EdgeTypeStr[INS] = "ins";
		EdgeTypeStr[DEL] = "del";
		EdgeTypeStr[SUBST] = "subst";
		EdgeTypeStr[JUMP] = "jump";
	}

	bool node_in_trie(int v) const {
		return v >= trie_first_node;
	}

	size_t trie_mem_bytes() const {
		return trie_edges * sizeof(E.front()) + trie_nodes * sizeof(V.front());
	}

	size_t total_mem_bytes() const {
		return E.size() * sizeof(E.front()) + V.size() * sizeof(V.front());
	}

	size_t total_mem_bytes_capacity() const {
		return E.capacity() * sizeof(E.front()) + V.capacity() * sizeof(V.front());
	}

	size_t reference_mem_bytes() const {
		return total_mem_bytes() - trie_mem_bytes();
	}

	int nodes() const {
		return V.size();
	}

	int edges() const {
		return E.size();
	}

  public:
	void init(int _n, int _m) {
		V.resize(_n, -1);  // 0 preserved for a supersource
		E.reserve(_m);
	}

	bool has_supersource() const {
		return V[0] != -1;
	}

	bool has_node(int u) const {
		return V[u] != -1;
	}

    int add_node() {
		V.push_back(-1);
		return V.size()-1;
	}

	void add_edge(int from, edge_t e) {
        E.push_back(e);
        V[from] = (int)E.size()-1;
	}

    void add_edge(int a, int b, char label, EdgeType type, int node_id=-1, int offset=-1) {
		LOG_FATAL_IF(!(a >= 0 && a < nodes())) << "edge with a=" << a << ", b=" << b << ", nodes=" << nodes();
		assert(a >= 0 && a < nodes());
		assert(b >= 0 && b < nodes());
		edge_t e(a, b, label, V[a], type, node_id, offset);
        E.push_back(e);
        V[a] = (int)E.size()-1;
    }

	void add_seq(int from, const std::string &seq, int to) {
		int prev=from;

		assert(seq.length() > 0);
		int curr = add_node();
		add_edge(prev, curr, seq[0], ORIG);
		prev = curr;

		for (std::string::size_type i=1; i<seq.size(); i++) {
			int curr = add_node();
			if (is_nucl(seq[i])) {
				add_edge(prev, curr, seq[i], ORIG);
			} else {
				assert(is_extended_nucl(seq[i]));
				std::string origs = extended2orignucls(seq[i]);
				for (char nucl: origs) {
					add_edge(prev, curr, nucl, ORIG);
				}
			}
			prev = curr;
		}
		add_edge(prev, to, EPS, ORIG);
	}

	void add_reverse_complement() {
#ifndef NDEBUG
		static bool done = false;
		assert(!done);
		assert(!has_supersource());
#endif

		// prepare the new nodes and edges
		int half_nodes = V.size();
		std::vector< std::pair<std::pair<int,int>, label_t> > new_edges;
		for (int from=0; from<(int)nodes(); from++) {
			for (int idx=V[from]; idx!=-1; idx=E[idx].next) {
				edge_t e = E[idx];
				new_edges.push_back(std::make_pair(std::make_pair(half_nodes + e.to, half_nodes + from), compl_nucl(e.label)));
			}
		}

		// add the new nodes and edges
		for (int i=0; i<half_nodes; i++)
			add_node();  // including the supersource which will stay unused
		for (const auto &e: new_edges)
			add_edge(e.first.first, e.first.second, e.second, astarix::ORIG);

#ifndef NDEBUG
		done = true;
#endif
	}

	void writeToStdout() const {
		printf("%d %d\n", (int)nodes(), (int)edges());
		for (int from=0; from<(int)nodes(); from++) {
			for (int idx=V[from]; idx!=-1; idx=E[idx].next) {
				edge_t e = E[idx];
				printf("%d %d %c %s\n", from, (int)e.to, (char)e.label, EdgeTypeStr[e.type]);
			}
		}
	}

	bool hasOutgoingEdges(int u) const {
		for (int idx=V[u]; idx!=-1; idx=E[idx].next) {
			if (E[idx].to != u)
				return true;
		}
		return false;
	}

	edge_t getOrigEdge(int u, int v) const {
		edge_t res;
		for (auto it=begin_orig_edges(u); it!=end_orig_edges(); ++it)
			if (it->to == v)
				res = *it;
		return res;
	}

	int numOutOrigEdges(int u, edge_t *e) const {
		int cnt=0;
		for (int idx=V[u]; idx!=-1; idx=E[idx].next)
			if (E[idx].type == ORIG) {
				cnt++;
				*e = E[idx];
			}
		return cnt;
	}

	//// ORIG EDGES ITERATOR (excl. edit edges)
	class orig_edge_iterator;
	orig_edge_iterator begin_orig_edges(int v) const { return orig_edge_iterator(this, v); }
	orig_edge_iterator end_orig_edges() const { return orig_edge_iterator(this, -1); }

	// Iterator of the original outgoing edges in the graph (excluding edit-edges).
	class orig_edge_iterator {
		const graph_t *g;
		int curr_edge_idx;

	  public:
		using value_type = edge_t;
		using reference = edge_t;
		using iterator_category = std::input_iterator_tag;
		using pointer = edge_t*;
		using difference_type = void;

		orig_edge_iterator(const graph_t *G, int _v)
			: g(G), curr_edge_idx(_v != -1 ? G->V[_v] : -1) {
		}

		const reference operator*() const { return g->E[curr_edge_idx]; }
		const pointer operator->() const { return (const pointer)&(g->E[curr_edge_idx]); }

		orig_edge_iterator& operator++() {  // preincrement
			curr_edge_idx = g->E[curr_edge_idx].next;
			return *this;
		}

		const orig_edge_iterator& operator++(int) { // postincrement
			const auto tmp = this;
			++(*this);
			return *tmp;
		}

		friend bool operator==(orig_edge_iterator const& lhs, orig_edge_iterator const& rhs) {
			return lhs.curr_edge_idx == rhs.curr_edge_idx;
		}

		friend bool operator!=(orig_edge_iterator const& lhs, orig_edge_iterator const& rhs) {
			return !(lhs == rhs);
		}
	};

	class all_matching_edge_iterator;

	all_matching_edge_iterator begin_all_edges(int v) const { return all_matching_edge_iterator(this, v, '!'); }
	all_matching_edge_iterator end_all_edges() const { return all_matching_edge_iterator(this, -1, '!'); }

	all_matching_edge_iterator begin_all_matching_edges(int v, label_t l) const { return all_matching_edge_iterator(this, v, l); }
	all_matching_edge_iterator end_all_matching_edges() const { return all_matching_edge_iterator(this, -1, '!'); }

	// Iterator of all outgoing edges in the graph (incl. edit-edges).
	class all_matching_edge_iterator {
		std::vector<edge_t> edit_edges;
		std::vector<edge_t>::const_iterator it;

	  public:
		using value_type = edge_t;
		using reference = edge_t;
		using iterator_category = std::input_iterator_tag;
		using pointer = edge_t*;
		using difference_type = void;

		all_matching_edge_iterator(const graph_t *G, int v, label_t l) {
			if (l != '!') {
				edit_edges.reserve(10);
				for (int idx=G->V[v]; idx!=-1; idx=G->E[idx].next) {
					const edge_t &orig_e = G->E[idx];
					if (orig_e.label == l)
						// match
						edit_edges.push_back(orig_e);
				}

				for (int idx=G->V[v]; idx!=-1; idx=G->E[idx].next) {
					const edge_t &orig_e = G->E[idx];
					if (orig_e.label != l)
						// substitution
						edit_edges.push_back(edge_t::from_cost(v, orig_e.to, l, SUBST));

					// deletions
					edit_edges.push_back(edge_t::from_cost(v, orig_e.to, EPS, DEL));
				}

				// insertions
				edit_edges.push_back(edge_t::from_cost(v, v, l, INS));
			}

			it = edit_edges.begin();
		}

		const reference operator*() const { return *it; }
		const pointer operator->() const { return (const pointer)&(*it); }

		all_matching_edge_iterator& operator++() {  // preincrement
			++it;
			return *this;
		}

		const edge_t operator++(int) { // postincrement
			edge_t tmp = *it;
			++it;
			return tmp;
		}

		friend bool operator==(all_matching_edge_iterator const& lhs, all_matching_edge_iterator const& rhs) {
			return (lhs.it == rhs.it) || (lhs.it == lhs.edit_edges.end() && rhs.it == rhs.edit_edges.end());
		}

		friend bool operator!=(all_matching_edge_iterator const& lhs, all_matching_edge_iterator const& rhs) {
			return !(lhs == rhs);
		}
	};
};

struct seq_t {
    std::string s;
    std::string comment;

	seq_t() {}
	seq_t(const std::string &_s, std::string _comment="")
		: s(_s), comment(_comment) {
		assert(are_all_nucls(s));
	}
};

struct read_t {
    std::string s; // 1-based
    std::string phreds;  // 1-based
    int len;  // s[1:len] has len letters
	//bool with_errors;
	std::string comment;
	std::string grnd_s;
	edge_path_t edge_path; //alignment

	cost_t mapping_quality;

	read_t() {}

	int size() const {
		return s.length() - 1;
	}

	read_t(const std::string &_s, std::string _phreds="", std::string _comment="", std::string _grnd_s="")
			: len(_s.length()), mapping_quality(-1) {
		assert(_s=="" || _s[0] != '@');
		assert(_grnd_s=="" || _grnd_s[0] != '+');

		comment = _comment;
		s = '@'	+ _s;
		grnd_s = '+' + _grnd_s;

		assert(are_all_nucls(_s));

		if (_phreds != "") {
			phreds = '@' + _phreds;
			assert(s.length() == phreds.length());
		}
	}
};

}

#endif
