#pragma once

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <plog/Log.h>
#include <plog/Appenders/ColorConsoleAppender.h>

#include "gfa2graph.h"
#include "graph.h"
#include "utils.h"

namespace astarix {

void _mkdir(const char *dir);
void assure_dir_exists(const char *dir_str);
std::string to_str(int argc, char **argv);

inline char uppercase(char c);
inline char lowercase(char c);

// Input
std::vector<seq_t> read_fasta(const std::string &fn);
void read_graph(graph_t *G, std::string graph_file, std::string output_dir);
bool read_query(std::ifstream &in, const std::string fn, read_t *r);

std::string spell(const edge_path_t &path);
std::string get_read_matching(const edge_path_t &path, const read_t &r);

// Output
void output_summary(const read_t &r, const EditCosts &costs, const edge_path_t &path, std::ostream &out);
void output_alignement(const graph_t &G, const EditCosts &costs, const read_t &r, const edge_path_t &path, std::ostream &out);
int output(const graph_t &G, const EditCosts &costs, const read_t &r, const edge_path_t &path, std::string output_file);

}
