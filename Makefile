

GPP=$(CXX)
CPPFLAGS=-Wall -DNDEBUG -Wextra -std=c++14 -O3 -Iext/plog/include/ -Iext/GraphAligner/ -Iext/concurrentqueue/ -Iext/parallel_hashmap/ -Wno-unused-parameter -Wno-missing-field-initializers

SRCDIR=src
EXTDIR=ext
DATADIR=data
TESTSDIR=tests

ODIR=obj
BINDIR=bin
TMPDIR=tmp

ASTARIXBIN=$(BINDIR)/astarix
LIBS= #-lm -lz 

_DEPS = $(SRCDIR)/argparse.h $(SRCDIR)/dijkstra.h $(SRCDIR)/astar-prefix.h $(SRCDIR)/astar-landmarks.h $(SRCDIR)/gfa2graph.h $(SRCDIR)/graph.h $(SRCDIR)/io.h $(SRCDIR)/align.h $(SRCDIR)/utils.h $(SRCDIR)/trie.h $(EXTDIR)/GraphAligner/GfaGraph.h
DEPS = $(patsubst %, %, $(_DEPS))

_OBJ = $(SRCDIR)/argparse.o $(SRCDIR)/astar-prefix.o $(SRCDIR)/gfa2graph.o $(SRCDIR)/align.o $(SRCDIR)/utils.o $(SRCDIR)/trie.o $(EXTDIR)/GraphAligner/GfaGraph.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

$(shell mkdir -p bin)
$(shell mkdir -p obj/src)
$(shell mkdir -p obj/ext/GraphAligner)

$(ASTARIXBIN): $(SRCDIR)/astarix.cpp $(OBJ)
	$(GPP) $< -o $@ $(OBJ) $(LINKFLAGS)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

test: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))

	# small
	$(ASTARIXBIN) align-optimal -a astar-prefix -t 8 -v 2 -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-prefix
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 1 -v 2 -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-landmarks --fixed_trie_depth 1
	$(ASTARIXBIN) align-optimal -a dijkstra -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/dijkstra-default
	python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head10000_linear/astar-prefix/alignments.tsv $(TMPDIR)/ecoli_head10000_linear/dijkstra-default/alignments.tsv
	python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head10000_linear/astar-landmarks/alignments.tsv $(TMPDIR)/ecoli_head10000_linear/dijkstra-default/alignments.tsv

bigtest:
	# 10000 reads
	$(ASTARIXBIN) align-optimal -t 8 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default
	$(ASTARIXBIN) align-optimal -a dijkstra -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/dijkstra-default
	python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head1000000_linear/astar-default/alignments.tsv $(TMPDIR)/ecoli_head1000000_linear/dijkstra-default/alignments.tsv

eval100: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 8 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default -v 0 --fixed_trie_depth 1
	$(ASTARIXBIN) align-optimal -a astar-prefix -t 8 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default -v 0 --fixed_trie_depth 1
                    
eval150: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 8 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina150.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default -v 0 --fixed_trie_depth 1
	$(ASTARIXBIN) align-optimal -a astar-prefix -t 8 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina150.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default -v 0 --fixed_trie_depth 1

eval250: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 8 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina250.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default -v 0 --fixed_trie_depth 1
	$(ASTARIXBIN) align-optimal -a astar-prefix -t 8 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina250.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default -v 0 --fixed_trie_depth 1

# pbsim --data-type CLR --depth 2 --model_qc ../pbsim_profiles/model_qc_clr ecoli.fasta
eval_long: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 8 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads/astar-default -v 0 --fixed_trie_depth 1

#pbsim --data-type CLR --depth 2 --model_qc ../pbsim_profiles/model_qc_clr --accuracy-mean 0.90 --length-min 300 --length-max 300 ecoli.fasta
eval_long300_90: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 8 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long300_90.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads/astar-default -v 0 --fixed_trie_depth 1

.PHONY: all clean

clean:
	rm -rf $(ODIR)/*
	rm -rf $(BINDIR)/*
	rm -rf $(TMPDIR)/*
