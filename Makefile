GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++14 -O3 -Iext/plog/include/ -Iext/GraphAligner/ -Iext/concurrentqueue/ -Iext/parallel_hashmap/ -Wno-unused-parameter -Wno-missing-field-initializers

DEBUG ?= 0
ifeq ($(DEBUG), 1)
    CPPFLAGS += -DDEBUG
	BINDIR=debug
	RUNFLAGS += -v 2
else
    CPPFLAGS += -DNDEBUG
	BINDIR=release
	RUNFLAGS += -v 0
endif

SRCDIR=src
EXTDIR=ext
DATADIR=data
TESTSDIR=tests

ODIR=$(BINDIR)/obj
TMPDIR=tmp

ASTARIXBIN=$(BINDIR)/astarix
LIBS= #-lm -lz 

_DEPS = $(SRCDIR)/argparse.h $(SRCDIR)/dijkstra.h $(SRCDIR)/astar-prefix.h $(SRCDIR)/astar-landmarks-exact.h $(SRCDIR)/astar-landmarks-approx.h $(SRCDIR)/gfa2graph.h $(SRCDIR)/graph.h $(SRCDIR)/io.h $(SRCDIR)/align.h $(SRCDIR)/utils.h $(SRCDIR)/trie.h $(EXTDIR)/GraphAligner/GfaGraph.h
DEPS = $(patsubst %, %, $(_DEPS))

_OBJ = $(SRCDIR)/argparse.o $(SRCDIR)/astar-prefix.o $(SRCDIR)/gfa2graph.o $(SRCDIR)/align.o $(SRCDIR)/utils.o $(SRCDIR)/trie.o $(EXTDIR)/GraphAligner/GfaGraph.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

$(shell mkdir -p $(BINDIR))
$(shell mkdir -p $(ODIR)/src)
$(shell mkdir -p $(ODIR)/ext/GraphAligner)

$(ASTARIXBIN): $(SRCDIR)/astarix.cpp $(OBJ)
	$(GPP) $< -o $@ $(OBJ) $(LINKFLAGS)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

test: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))

	# small
	$(ASTARIXBIN) align-optimal -a astar-landmarks-exact -t 1 $(RUNFLAGS) -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-landmarks-exact --fixed_trie_depth 1 --astar_max_waymark_errors 0
	$(ASTARIXBIN) align-optimal -a astar-landmarks       -t 1 $(RUNFLAGS) -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-landmarks-approx --fixed_trie_depth 1 --astar_max_waymark_errors 0
	$(ASTARIXBIN) align-optimal -a astar-landmarks       -t 1 $(RUNFLAGS) -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-landmarks-approx --fixed_trie_depth 1 --astar_max_waymark_errors 1
	#$(ASTARIXBIN) align-optimal -a astar-prefix -t 6 $(RUNFLAGS) -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-prefix
	#$(ASTARIXBIN) align-optimal -a dijkstra -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/dijkstra-default

	#python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head10000_linear/astar-prefix/alignments.tsv $(TMPDIR)/ecoli_head10000_linear/dijkstra-default/alignments.tsv
	#python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head10000_linear/astar-landmarks-exact/alignments.tsv $(TMPDIR)/ecoli_head10000_linear/dijkstra-default/alignments.tsv
	#python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head10000_linear/astar-landmarks-approx/alignments.tsv $(TMPDIR)/ecoli_head10000_linear/dijkstra-default/alignments.tsv

bigtest:
	# 10000 reads
	$(ASTARIXBIN) align-optimal -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default
	$(ASTARIXBIN) align-optimal -a dijkstra -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/dijkstra-default
	python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head1000000_linear/astar-default/alignments.tsv $(TMPDIR)/ecoli_head1000000_linear/dijkstra-default/alignments.tsv

eval100: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default $(RUNFLAGS) --fixed_trie_depth 1
	$(ASTARIXBIN) align-optimal -a astar-prefix -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default $(RUNFLAGS) --fixed_trie_depth 1
                    
eval150: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina150.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default $(RUNFLAGS) --fixed_trie_depth 1
	$(ASTARIXBIN) align-optimal -a astar-prefix -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina150.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default $(RUNFLAGS) --fixed_trie_depth 1

eval250: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina250.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default $(RUNFLAGS) --fixed_trie_depth 1
	$(ASTARIXBIN) align-optimal -a astar-prefix -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina250.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default $(RUNFLAGS) --fixed_trie_depth 1

# pbsim --data-type CLR --depth 2 --model_qc ../pbsim_profiles/model_qc_clr ecoli.fasta
eval_long: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	#$(ASTARIXBIN) align-optimal -a astar-landmarks-exact -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_max_waymark_errors 0
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_max_waymark_errors 2

#pbsim --data-type CLR --depth 2 --model_qc ../pbsim_profiles/model_qc_clr --accuracy-mean 0.90 --length-min 300 --length-max 300 ecoli.fasta
eval_long300_90: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long300_90.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads/astar-waymark_errors_2 $(RUNFLAGS) --fixed_trie_depth 1 --astar_max_waymark_errors 3 -G 1 -S 1 --astar_waymark_len 15 --tree_depth 13

eval_long300_05: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long300_05.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads/astar-waymark_errors_2 $(RUNFLAGS) --fixed_trie_depth 1 --astar_max_waymark_errors 2 -G 1 -S 1 --astar_waymark_len 15 --tree_depth 13

eval_manual: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-landmarks -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long_manual.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads_manual/astar-waymark_errors_2 $(RUNFLAGS) --fixed_trie_depth 1 --astar_max_waymark_errors 2 -G 1 -S 1 

.PHONY: all clean

clean:
	#rm -rf $(ODIR)/*
	rm -rf debug/
	rm -rf release/
	rm -rf $(TMPDIR)/
