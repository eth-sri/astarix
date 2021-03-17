#include "graph.h"

namespace astarix {

std::ostream& operator<<(std::ostream& os, const state_t &st) {
    os << "(i=" << st.i << ", v=" << st.v << ", cost=" << st.cost << ")";
    return os;
}

}