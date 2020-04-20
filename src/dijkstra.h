#pragma once

#include <map>
#include <string>
#include <vector>

#include "graph.h"
#include "utils.h"

namespace astarix {

class DijkstraDummy: public AStarHeuristic {

  public:

    cost_t h(const state_t &st) const {
        return 0;
    }

    void before_every_alignment(const read_t *r) {
        // Empty.
    }

    void after_every_alignment() {
        // Empty.
    }

    void print_params(std::ostream &out) const {
        // Empty.
    }

    void print_stats(std::ostream &out) const {
        // Empty.
    }
};

}
