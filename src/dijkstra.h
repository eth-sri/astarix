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
};

}
