GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -DNDEBUG -std=c++14 -O3 -Iext/plog/include/ -Iext/GraphAligner/ -Iext/concurrentqueue/ -Iext/parallel_hashmap/ -Wno-unused-parameter

SRCDIR=src
EXTDIR=ext
DATADIR=data
TESTSDIR=tests

ODIR=obj
BINDIR=bin
TMPDIR=tmp

ASTARIXBIN=$(BINDIR)/astarix
LIBS= #-lm -lz 

_DEPS = $(SRCDIR)/argparse.h $(SRCDIR)/astar.h $(SRCDIR)/gfa2graph.h $(SRCDIR)/graph.h $(SRCDIR)/io.h $(SRCDIR)/align.h $(SRCDIR)/utils.h $(SRCDIR)/trie.h $(EXTDIR)/GraphAligner/GfaGraph.h
DEPS = $(patsubst %, %, $(_DEPS))

_OBJ = $(SRCDIR)/argparse.o $(SRCDIR)/astar.o $(SRCDIR)/gfa2graph.o $(SRCDIR)/align.o $(SRCDIR)/trie.o $(EXTDIR)/GraphAligner/GfaGraph.o
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
	$(ASTARIXBIN) align-optimal -t 8 -v 2 -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-default
	$(ASTARIXBIN) align-optimal -a dijkstra -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/dijkstra-default
	python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head10000_linear/astar-default/alignments.tsv $(TMPDIR)/ecoli_head10000_linear/dijkstra-default/alignments.tsv

bigtest:
	# 10000 reads
	$(ASTARIXBIN) align-optimal -t 8 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default
	$(ASTARIXBIN) align-optimal -a dijkstra -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/dijkstra-default
	python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head1000000_linear/astar-default/alignments.tsv $(TMPDIR)/ecoli_head1000000_linear/dijkstra-default/alignments.tsv

eval: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -t 8 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default
                    
.PHONY: all clean

clean:
	rm -rf $(ODIR)/*
	rm -rf $(BINDIR)/*
	rm -rf $(TMPDIR)/*
