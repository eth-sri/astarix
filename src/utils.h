#ifndef ASTARIX_UTILS_H
#define ASTARIX_UTILS_H

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <ctime>
#include <fstream>
#include <ios>
#include <iostream>
#include <limits.h>
#include <string>
#include <unistd.h>
#include <vector>
#include <queue>

namespace astarix {

#define eps (1e-6)
#define MIN(a,b) (((a)<(b))?(a):(b))
#define EQ(a, b) ((a<b) ? (b-a) < eps : (a-b) < eps)

typedef char cost_t;
typedef char label_t;

class state_t;
typedef std::pair<cost_t, state_t> 										score_state_t;
typedef std::priority_queue<score_state_t, std::vector<score_state_t>, std::greater<score_state_t>>
																			queue_t;
const char nucls[] = "ACGT";
const std::string extended_nucls = "RYKMSWBDHVN";
const char EPS   = 'e';

//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0
// https://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-runtime-using-c

static void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entrees in stat that we don't care about
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

static double b2gb(size_t bytes) {
	return bytes / 1024.0 / 1024.0 / 1024.0;
}

class Timer {
	clock_t start_time;
	clock_t accum_time;
	bool running;

  public:
	Timer() : accum_time(0), running(false) {}

	void clear() {
		accum_time = 0;
		running = false;
	}

	void start() {
		assert(!running);
		running = true;
		start_time = std::clock();
	}

	void stop() {
		accum_time += std::clock() - start_time;
		assert(running);
		running = false;
	}

	double get_sec() const {
		assert(!running);
		return (double)accum_time / (double)CLOCKS_PER_SEC;
	}

	Timer& operator+=(const Timer &b) {
		accum_time += b.accum_time;
		return *this;
	}
};

class MemoryMeasurer {
	double start_mem;
	double accum_mem;
	bool running;

  public:
	MemoryMeasurer() : accum_mem(0.0), running(false) {}

	void clear() {
		accum_mem = 0.0;
		running = false;
	}

	static double get_mem_gb() {
		double vm, rss;
		process_mem_usage(vm, rss);
		return rss / (1024.0 * 1024.0);  // KB to GB
		//return vm / (1024.0 * 1024.0);  // KB to GB
	}

	void start() {
		assert(!running);
		running = true;
		start_mem = get_mem_gb();
	}

	void stop() {
		accum_mem += get_mem_gb() - start_mem;
		assert(running);
		running = false;
	}

	double get_gb() const {
		assert(!running);
		return accum_mem;
	}
};

class Counter {
	typedef int T;
	T cnt;

  public:
	Counter() : cnt(T(0)) {}

	void clear() {
		cnt = 0;
	}

	void inc() {
		++cnt;
	}

	void inc(T a) {
		cnt += a;
	}

	T get() const {
		return cnt;
	}
};
static std::string bool2str(bool x) {
	return x ? "true" : "false";
}

enum EdgeType : char {
	ORIG, INS, DEL, SUBST, JUMP, EdgeType_after_type
};

typedef int nodesz;

struct edge_t {
    nodesz to;    // 4 bytes, the end point of the edge
	nodesz next;  // 4 bytes, E[next] -- prev added outgoing edge from the same source
    label_t label; // 1 byte automata label
	EdgeType type;  // 1 byte

    edge_t() : to(-1), next(-1), label(EPS), type(ORIG) {}
	edge_t(int _from, int _to, label_t _label, int _next, EdgeType _type=ORIG, int _node_id=-1, int _offset=-1)
		: to(_to), next(_next), label(_label), type(_type) {}

	static edge_t from_cost(int _from, int _to, label_t _label, EdgeType _type) {
		edge_t e;
		e.to = _to;
		e.label = _label;
		e.type = _type;
		return e;
	}
};

static inline const char* edgeType2str(EdgeType type) {
	switch(type) {
		case ORIG:  return "ORIG";
		case INS:   return "ins";
		case DEL:   return "del";
		case SUBST: return "subst";
		case JUMP:  return "jump";
		case EdgeType_after_type: assert(false);
	}
	return "?";
}

class EditCosts {
  public:
	cost_t match, subst, ins, del;
	EditCosts() : match(-1), subst(-1), ins(-1), del(-1) {}
	EditCosts(cost_t _match, cost_t _subst, cost_t _ins, cost_t _del)
		: match(_match), subst(_subst), ins(_ins), del(_del) {
	}
		
	cost_t edge2score(const edge_t &e) const {
		switch (e.type) {
			case ORIG:
			case JUMP:
				return match;
			case SUBST:
				return subst;
			case INS:
				return ins;
			case DEL:
				return del;
			case EdgeType_after_type:
				assert(false);
		}
		std::cerr << "Edge type: " << e.type << std::endl;
		assert(false && "No such edge type");
		return -10000.0;
	}
};

static inline char compl_nucl(char c) {
	if (c == 'A') return 'T';
	if (c == 'C') return 'G';
	if (c == 'G') return 'C';
	if (c == 'T') return 'A';
	std::cerr << "Bad nucleotide '" << c << "'" << std::endl;
	assert(false);
	return '!';
}

static inline int nucl2num(char c) {
	if (c == 'A') return 0;
	if (c == 'C') return 1;
	if (c == 'G') return 2;
	if (c == 'T') return 3;
	assert(false);
	return '!';
}

static inline bool is_nucl(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N';
}

static std::string to_lower(std::string s) {
	for (char &c: s)
		if (c >= 'A' && c <= 'Z')
			c -= 'A'-'a';
	return s;
}

static inline bool is_letter(char c) {
	return is_nucl(c) || c == EPS;
}

#ifndef NDEBUG
static bool is_extended_nucl(char nucl) {
	return extended_nucls.find(nucl) != std::string::npos;
}

static bool are_all_nucls(const std::string &s) {
	for (char c: s) {
		if (!is_nucl(toupper(c))) {
			std::cerr << "Not a nucleotide: [" << c << "]" << std::endl;
			return false;
		}
	}
	return true;
}

static void write(std::ostream& os, int from, const edge_t &e) {
	os << from << " " << e.to << " " << e.label << " " << edgeType2str(e.type); 
}
#endif

static std::string extended2orignucls(char nucl) {
	if (nucl == 'R') return "AG";
	if (nucl == 'Y') return "CT";
	if (nucl == 'K') return "GT";
	if (nucl == 'M') return "AC";
	if (nucl == 'S') return "CG";
	if (nucl == 'W') return "AT";
	if (nucl == 'B') return "CGT";
	if (nucl == 'D') return "AGT";
	if (nucl == 'H') return "ACT";
	if (nucl == 'V') return "ACG";
	if (nucl == 'N') return "ACGT";
	assert(is_nucl(nucl));
	return "!";
}

static inline double sample() {
	return 1.0 * rand() / INT_MAX;
}

// https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
static bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

}

#endif
