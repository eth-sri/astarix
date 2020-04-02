#pragma once

#include <string>

#include <GfaGraph.h>
#include "graph.h"

GfaGraph load_gfa(const std::string &gfa_filename);
void gfa2graph(GfaGraph &gfa, astarix::graph_t *G);
