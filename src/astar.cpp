#include "astar.h"

namespace astarix {

bool AStar::is_linear(int u, int rem_len, std::string *pref, int *boundary_node) {
	if (rem_len == 0) {
		(*boundary_node) = u;
		return true;
	}

	int v=-1;
	char c;
	for (auto it=G.begin_orig_edges(u); it!=G.end_orig_edges(); ++it) {
		const edge_t &e = *it;
		if (v != -1)
			return false;
		v = e.to;
		c = e.label;
	}

	if (v == -1)       									// zero neighbours
		return false;

	(*pref) += c;
	return is_linear(v, rem_len-1, pref, boundary_node);
}

int AStar::precompute_A_star_prefix() {
	LOG_INFO << "A* precomputation...";
	assert(_star.empty());

	LOG_INFO << "Using A* prefix len " << max_prefix_len << " and A* max cost " << max_prefix_cost;
	_vertex2class.resize(G.nodes());

	hash_precomp();  									// computed kMaxStrHash
	assert(kMaxStrHash != -1);

	std::vector<int> strhash2class(kMaxStrHash, -1);  	// internal hashing: a hash of a prefix of a linear part of the graph to a representative vertex with a computed future

	classes = 0;

	// nodes compressed to classes of equivalence
	for (int u=0; u<G.nodes(); u++) {
		std::string pref;
		int boundary_node=-1;
		if (compress_vertices && is_linear(u, max_prefix_len, &pref, &boundary_node)) {
			++compressable_vertices;
			auto h = hash_str(pref);
			assert((size_t)h < strhash2class.size());
			int &cl = strhash2class[h];   			// representative vertex
			if (cl == -1) {  						// if this is the first vertex with seen from this class
				cl = classes++;
				_class2repr.push_back(u);
				_class2boundary.push_back(boundary_node);
			}
			_vertex2class[u] = cl;   				// to be used for quering
		} else {
			int cl = classes++;
			_class2repr.push_back(u);
			_class2boundary.push_back(-1);
			_vertex2class[u] = cl;   
		}
	}

	LOG_INFO << G.nodes() << " vertices compressed to " << classes
		<< " A* equivalence classes (" << 1.0*G.nodes()/classes << "x memory saving)";
	LOG_INFO << "Hash table size: " << size_t(classes)*kMaxStrHash << " = " << classes << " * " << kMaxStrHash;

	LOG_INFO << "Bucket count before resizing: " << _star.bucket_count();
	//_star.resize(150 * 1000 * 1000);
	//LOG_INFO << "Bucket count after initial resizing: " << _star.bucket_count();
	assert(_class2repr.size() == (size_t)classes);
	LOG_INFO << "Prefix+Vertex hash table size = " << _star.size();

	LOG_INFO << "Precomputation finished.";

	int precomputed_elements = (int)_star.size();
	return precomputed_elements;
}

void AStar::compute_astar_cost_from_vertex_and_prefix(cost_t &res, int u, const std::string &prefix,
													   int boundary_node, int i, cost_t prev_cost) {
	if ((size_t)i >= prefix.size() || u == boundary_node) {
		if (prev_cost < res) {
			res = prev_cost;
		}
		return;
	}

	if (prev_cost >= res)
		return;

	for (auto it=G.begin_all_matching_edges(u, prefix[i]); it!=G.end_all_matching_edges(); ++it) {
		const edge_t &e = *it;
		int next_i = e.label == EPS ? i : i+1;
		compute_astar_cost_from_vertex_and_prefix(res, e.to, prefix, boundary_node, next_i, prev_cost+costs.edge2score(e));
	}
}

cost_t AStar::lazy_star_value(unsigned h, int repr, int boundary_node, const std::string &prefix) {
	LOG_DEBUG << "Lazy A* query for h=" << h << ", repr=" << repr << ", boundary_node=" << boundary_node << ", prefix=" << prefix;

	++_cache_trees;
	auto it = _star.find(h);
	if (it == _star.end()) {
		++_cache_misses;
		cost_t res=0.0;
		res = max_prefix_cost;
		Timer curr_time;
		curr_time.start();
			compute_astar_cost_from_vertex_and_prefix(res, repr, prefix, boundary_node);
		curr_time.stop();
		astar_time += curr_time;
		++_entries;
		_star[h] = res;
		return res;
	}

	return it->second;
}

cost_t AStar::astar_from_pos(int v, const std::string &prefix) {
	LOG_DEBUG << "v=" << v << ", prefix=" << prefix;
	assert(v < _vertex2class.size());
	int cl = _vertex2class[v];
	auto h = hash(prefix, cl);
	assert(cl < _class2repr.size());
	int repr = _class2repr[cl];
	assert(cl < _class2boundary.size());
	int boundary_node = _class2boundary[cl];
	return lazy_star_value(h, repr, boundary_node, prefix);
}

cost_t AStar::h(int v, const std::string &prefix) {
	LOG_FATAL_IF(prefix.length() > (size_t)max_prefix_len)
		<< "The prefix " << prefix << " with length " << prefix.length() << " should be shorter than " << max_prefix_len;
	assert(prefix.length() <= (size_t)max_prefix_len);
	return astar_from_pos(v, prefix);
}

}
