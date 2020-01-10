#ifndef ASTARIX_IO_H
#define ASTARIX_IO_H

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <dirent.h>
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

// http://nion.modprobe.de/blog/archives/357-Recursive-directory-creation.html
static void _mkdir(const char *dir) {
	char tmp[256];
	char *p = NULL;
	size_t len;

	snprintf(tmp, sizeof(tmp),"%s",dir);
	len = strlen(tmp);
	if(tmp[len - 1] == '/')
		tmp[len - 1] = 0;
	for(p = tmp + 1; *p; p++)
		if(*p == '/') {
			*p = 0;
			mkdir(tmp, S_IRWXU);
			*p = '/';
		}
	mkdir(tmp, S_IRWXU);
}

static void assure_dir_exists(char *dir_str) {
	DIR* dir = opendir(dir_str);
	if (dir) {  // Directory exists
		LOG_INFO << "Using dir " << dir_str;
		closedir(dir);
	}
	else if (ENOENT == errno) {  // Directory does not exist
		_mkdir(dir_str);
	}
	else {  // opendir() failed for some other reason
		LOG_ERROR << "Error opening dir " << dir_str;
		exit(0);
	}
}

static std::string to_str(int argc, char **argv) {
	std::string s;
	for(int i=0; i<argc; i++)
		s += argv[i] + std::string(" ");
	return s;
}

static std::vector<seq_t> read_fasta(const std::string &fn) {
	std::vector<seq_t> res;
	std::string comment;
	std::ifstream in(fn, std::ifstream::in);
	if(in.fail()) {
		LOG_FATAL << "Problem opening fasta file " << fn;
		assert(false);
	}

	getline(in, comment);
	while(true) {
		assert(comment != "");
		std::string line, s;
		while(getline(in, line), line.size() && line[0]!='>') {
			assert(line != "");
			s += line;
		}

		if (s.empty())
			break;

		std::transform(s.begin(), s.end(), s.begin(), ::toupper);
		seq_t seq(s, comment.substr(1));
		res.push_back(seq);
		comment = s;
		LOG_INFO << "Fasta record with length " << s.size() << " read from the file " << fn;
	}
	return res;
}

static void read_graph(graph_t *G, std::string graph_file, std::string output_dir) {
	LOG_INFO << "Reading graph " << graph_file << "...";

	if (hasEnding(to_lower(graph_file), ".fa") || hasEnding(to_lower(graph_file), ".fasta")) {
		G->add_node();  // for supersource

		std::vector<astarix::seq_t> fastas = astarix::read_fasta(graph_file);
		for (const auto &fasta: fastas) {
			int source=G->add_node();
			int sink=G->add_node();
			G->add_seq(source, fasta.s, sink);
		}
	} else if (hasEnding(to_lower(graph_file), ".gfa")) {
		LOG_INFO << "[GFA format]";
		GfaGraph gfa = load_gfa(graph_file);
		gfa2graph(gfa, G);
		G->add_reverse_complement();
	} else {
		LOG_INFO << "[unknown format]";
		assert(false);
	}

	assert(!G->has_supersource());
}

static bool read_query(std::ifstream &in, const std::string fn, read_t *r) {
	std::string s, comment, grnd_s;
	if (in.is_open()) {
		if (!getline(in, comment))
			return false;
		getline(in, s);
		std::transform(s.begin(), s.end(), s.begin(), ::toupper);
		std::for_each(s.begin(), s.end(), [](char c) { assert(is_nucl(c)); });
	} else {
	    LOG_ERROR << "Cannot read query file " << fn;
	    assert(false);
	}

	std::string phreds;

	if (hasEnding(fn, ".fastq") || hasEnding(fn, ".fq")) {
		getline(in, grnd_s);
		getline(in, phreds);
		LOG_ERROR_IF(grnd_s[0] != '+') << "grnd_s = " << grnd_s << " does not start with +";
		assert(grnd_s[0] == '+');
		assert(!comment.empty() && comment[0] == '@');
	} else {
		LOG_ERROR_IF(!hasEnding(fn, ".fasta") && !hasEnding(fn, ".fa")) << fn << " should be a .fasta/.fa or .fastq/.fq file.";
		grnd_s = "+";
	}
	*r = read_t(s, phreds, comment.substr(1), grnd_s.substr(1));

	assert(r->len > 0);

	return true;
}

static inline char uppercase(char c) { return c<='Z' ? c : (c-('a'-'A')); }
static inline char lowercase(char c) { return c>='a' ? c : (c+('a'-'A')); }

static std::string spell(const edge_path_t &path) {
	std::string ans;
	for (const auto &e: path) {
		if (e.label != EPS) {
			char label = (e.type == ORIG || e.type == JUMP) ? uppercase(e.label) : lowercase(e.label);
			ans.push_back(label);
		}
	}
	return ans;
}

static std::string get_read_matching(const edge_path_t &path, const read_t &r) {
	int last_read_i=1;
	std::string read_match;
	for (edge_t curr: path) {
		if (curr.label != EPS) {
			assert(last_read_i > 0 && last_read_i <= r.len);
			read_match += (curr.label == r.s[last_read_i]) ? '-' : r.s[last_read_i];
			last_read_i++;
		} else {
			read_match += '.';
		}
	}
	return read_match;
}

static void output_summary(const read_t &r, const EditCosts &costs, const edge_path_t &path, std::ostream &out) {
	out << " -- Summary -- " << std::endl;

	std::map<EdgeType, cost_t> type2scoresum;
	std::map<EdgeType, int> type2cnt;

	for (int i=0; i<EdgeType_after_type; i++) {
		type2scoresum[static_cast<EdgeType>(i)] = 0.0;
		type2cnt[static_cast<EdgeType>(i)] = 0;
	}

	for (edge_t curr: path) {
		type2scoresum[curr.type] += costs.edge2score(curr);
		type2cnt[curr.type]++;
	}

	for(auto iter: type2scoresum) {
		out << "    " << type2cnt[iter.first] << " " << edgeType2str(iter.first) << " edges sum to a cost of " << iter.second << std::endl;
	}
	std::string m = get_read_matching(path, r);
	out << "  " << m.size() - count(m.begin(), m.end(), '-') - count(m.begin(), m.end(), '.') << std::endl;
	LOG_DEBUG << "path: " << spell(path);
	LOG_DEBUG << "read: " << r.s.substr(1);

	out << std::endl;
}

static void output_alignement(const graph_t &G, const EditCosts &costs, const read_t &r, const edge_path_t &path, std::ostream &out) {
	std::string grnd_s, r_phred, read, align_path, orig_path, t, edge_phred, read_match;
	int read_idx=1;

	std::string read_phreds = r.phreds != "" ? r.phreds : std::string(r.s.length(), '?');

	for (edge_t curr: path) {
		char symb, read_phred, grnd_c;
		if (curr.label != EPS) {
			symb = r.s[read_idx];
			read_phred = read_phreds[read_idx];
			grnd_c = r.grnd_s[read_idx];
			read_idx++;
		} else {
			symb = '.';
			read_phred = '.';
			grnd_c = '.';
		}

		grnd_s += grnd_c;
		r_phred += read_phred;
		read += symb;
		align_path += curr.label;
		orig_path += G.getOrigEdge(curr.to, curr.to).label;  // TODO: should be curr.from
		t += edgeType2str(curr.type)[0] != 'O' ? edgeType2str(curr.type)[0] : '.';
	}
	read_match = get_read_matching(path, r);

	assert(read_match.size() == r_phred.size());
	assert(read.size() == r_phred.size());
	assert(align_path.size() == r_phred.size());
	assert(orig_path.size() == r_phred.size());
	assert(t.size() == r_phred.size());

	out << " -- Alignment -- " << std::endl;
	out << "read qual: " << r_phred << std::endl;
	out << "read erro: " << read_match << std::endl;
	out << "orig read: " << read << std::endl;
	out << "alig path: " << align_path << std::endl;
	out << "orig path: " << orig_path << std::endl;
	out << "grnd seq?: " << grnd_s << std::endl;
	out << "edge type: " << t << std::endl;
	out << std::endl;
}

static int output(const graph_t &G, const EditCosts &costs, const read_t &r, const edge_path_t &path, std::string output_file) {
	std::ofstream out(output_file);
	if (!out) {
		LOG_ERROR << "Cannot open for writing the file " << output_file;
		return 1;
	}

	if (path.empty()) {
		LOG_ERROR << "Empty alignment path.\n";
		return -1;
	}

	// to out
	output_summary(r, costs, path, out);
	output_alignement(G, costs, r, path, out);

	return 0;
}

}

#endif
