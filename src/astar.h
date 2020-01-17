#ifndef ASTARIX_ASTAR_H
#define ASTARIX_ASTAR_H

#include <string>
#include <vector>
#include "phmap.h"  // parallel_hashmap

#include "utils.h"
#include "io.h"

namespace astarix {

class AStar {
  private:
	// General params
	const graph_t &G;
	const EditCosts &costs;

	int max_prefix_len;
	cost_t max_prefix_cost;
	int kMaxStrHash;  								// calculated here, in hash_precomp()
	bool compress_vertices;
	bool lazy;

	// Main data struct used for memoization
	phmap::parallel_flat_hash_map<unsigned, cost_t,
							phmap::container_internal::hash_default_hash<unsigned>, \
                            phmap::container_internal::hash_default_eq<unsigned>, \
                            std::allocator<std::pair<const unsigned, cost_t>>, 4, std::mutex> _star;  // with a mutex

	// Auxiliary structs
	std::vector<unsigned> _prev_group_sum;  		// string length -> number of strings with strictly lower length
	std::vector<int> _vertex2class; 				// vertex to a representative equivalent vertex for which future is calculated
	std::vector<int> _class2repr;  					// used for `decompresion'
	std::vector<int> _class2boundary;
	unsigned _nucl_num[256];
	int _cache_trees, _cache_misses;

	int classes;  									// number of equivalence classes

	// Currently calculated
	int compressable_vertices;

	Timer astar_time;
	std::atomic<int> _entries;

  public:
	AStar(const graph_t &_G, const EditCosts &_costs, int _max_prefix_len, int _max_prefix_cost, bool _compress_vertices)
		: G(_G),
		  costs(_costs),
		  max_prefix_len(_max_prefix_len),
		  max_prefix_cost(_max_prefix_cost),
		  compress_vertices(_compress_vertices),
		  lazy(true),
		  _cache_trees(0), _cache_misses(0),
		  classes(0),
		  compressable_vertices(0),
		  _entries(0)
	{
		LOG_INFO << "A* Prefix class constructed with:";
		LOG_INFO << "  max_prefix_len    = " << max_prefix_len;
		LOG_INFO << "  max_prefix_cost  = " << max_prefix_cost;
		LOG_INFO << "  compress_vertices = " << (compress_vertices ? "true" : "false");
		LOG_INFO << "  lazy              = " << (lazy ? "true" : "false");

		// used for hashing
		_nucl_num['a'] = _nucl_num['A'] = 0;
		_nucl_num['c'] = _nucl_num['C'] = 1;
		_nucl_num['g'] = _nucl_num['G'] = 2;
		_nucl_num['t'] = _nucl_num['T'] = 3;

		precompute_A_star_prefix();
	}

	// heuristic h(node, upcoming string)
	cost_t h(int v, const std::string &prefix);

	void reset_astar_time() {
		astar_time.clear();
	}

	int get_max_prefix_len() const {
		return max_prefix_len;
	}

	size_t equiv_classes_mem_bytes() const {
		return _vertex2class.size() * sizeof(_vertex2class.front()) +
			_class2repr.size() * sizeof(_class2repr) +
			_class2boundary.size() * sizeof(_class2boundary);
	}

	double get_astar_time() const {
		return astar_time.get_sec();
	}

	int get_compressable_vertices() const {
		return compressable_vertices;
	}

	cost_t get_max_prefix_cost() const {
		return max_prefix_cost;
	}

	int get_cache_trees() const {
		return _cache_trees;
	}

	int get_cache_misses() const {
		return _cache_misses;
	}

	double table_entrees() const {
		return _star.size();
	}

	size_t entries() {
		return _entries;
	}

	size_t table_mem_bytes_lower() const {
		return _star.size() * (sizeof(unsigned) + sizeof(cost_t) + 1) / _star.load_factor();
	}

	size_t table_mem_bytes_upper() const {
		size_t add_size = 0.03 *_star.size() * (sizeof(unsigned) + sizeof(cost_t) + 1) / 0.4375;
		return table_mem_bytes_lower() + add_size;
	}

  private:
	// returns the total number of precomputed elements.
	// space for paths: O(NlogL*?), time for paths: O(MlogL*?)
	// query will be for O(logL)
    int precompute_A_star_prefix();

	// returns true if there is a unique ORIG path from u with length rem_len; postcond: pref is the spelling of this path
	// returns false otherwise
	bool is_linear(int u, int rem_len, std::string *pref, int *boundary_node);

	// computes the minimum cost for matching the remining read (prefix[i:]) from u
	void compute_astar_cost_from_vertex_and_prefix(cost_t &res, int u, const std::string &prefix, int boundary_node,
			int i=0, cost_t prev_cost=0.0);

	// wrapper around compute_astar_cost_from_vertex_and_prefix dealing with memoization
	cost_t lazy_star_value(unsigned h, int repr, int boundary_node, const std::string &prefix);

	// translatex (node, prefix) to (hash(eq_class_representative_node(node)), prefix)
	cost_t astar_from_pos(int v, const std::string &prefix);

	void hash_precomp() {
		int four_power=1;
		_prev_group_sum.clear();
		_prev_group_sum.reserve(max_prefix_len+1);
		_prev_group_sum.push_back(0);
		for (int i=1; i<=max_prefix_len; i++) {
			_prev_group_sum.push_back(_prev_group_sum.back() + four_power);
			four_power <<= 2;
			assert(four_power >= 0);
		}	
		kMaxStrHash = _prev_group_sum[max_prefix_len] + four_power;  // accounting for the last group as well
		LOG_INFO << "kMaxStrHash = " << kMaxStrHash;
	}

	// returns [0; 4^0+4^1+...+4^max_prefix_len) = [0; kMaxStrHash)
	unsigned hash_str(const std::string &s) {
		unsigned h = 0;
		for (char c: s) {
			assert(is_nucl(c));
			h = (h<<2) + _nucl_num[(unsigned)c];
		}
		h += _prev_group_sum[s.length()];
		//LOG_INFO << "h(" << s << ") = " << h;
		return h;
	}

	// returns [0; N*2^(2*K+1))
	unsigned hash(const std::string &s, unsigned cl) {
		unsigned h = hash_str(s) + cl*kMaxStrHash;
		return h;
	}
};

}

#endif
