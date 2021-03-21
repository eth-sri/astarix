#pragma once

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

#define _unused(x) ((void)(x))   // to not have warnings

#define eps (1e-6)
#define MIN(a,b) (((a)<(b))?(a):(b))
#define EQ(a, b) ((a<b) ? (b-a) < eps : (a-b) < eps)

typedef int cost_t;
typedef char label_t;

class state_t;
typedef std::pair<cost_t, state_t>                                      score_state_t;
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

void process_mem_usage(double& vm_usage, double& resident_set);
double b2gb(size_t bytes);  // bytes to gb

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

template<typename T = int>
class Counter {
    T cnt;

  public:
    Counter() : cnt(T(0)) {}

    void clear() {
        cnt = T(0);
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

    void set(T a) {
        cnt = a;
    }

    T operator++() {
        return ++cnt;
    }

    Counter& operator+=(const Counter &b) {
        cnt += b.cnt;
        return *this;
    }

//    friend std::ostream& operator<<(std::ostream& os, const Counter& c);
    friend std::ostream& operator<<(std::ostream& os, const Counter& c) {
        os << c.cnt;
        return os;
    }
};

inline std::string bool2str(bool x) {
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

class EditCosts {
  public:
    cost_t match, subst, ins, del;
    EditCosts() : match(-1), subst(-1), ins(-1), del(-1) {}
    EditCosts(cost_t _match, cost_t _subst, cost_t _ins, cost_t _del)
        : match(_match), subst(_subst), ins(_ins), del(_del) {
    }

    cost_t get_min_mismatch_cost() const {
        return std::min(subst, std::min(ins, del));
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

const char* edgeType2str(EdgeType type);

char compl_nucl(char c);
int nucl2num(char c);
bool is_nucl(char c);
std::string to_lower(std::string s);
bool is_letter(char c);

// for debugging
bool is_extended_nucl(char nucl);
bool are_all_nucls(const std::string &s);
void write(std::ostream& os, int from, const edge_t &e);

// from a letter denoting a subset to a list of nucleotides
std::string extended2orignucls(char nucl);

// in [0,1]
double sample();

// https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
bool hasEnding(const std::string &fullString, const std::string &ending);

}
