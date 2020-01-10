GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++14 -O3 -g -Iext/plog/include/ -Iext/GraphAligner/ -Iext/ -Wno-unused-parameter

SRCDIR=src
EXTDIR=ext
DATADIR=data

ODIR=obj
BINDIR=bin
TMPDIR=tmp

ASTARIX_BIN=$(BINDIR)/astarix
LIBS=-lm -lz 

_DEPS = $(SRCDIR)/argparse.h $(SRCDIR)/astar.h $(SRCDIR)/gfa2graph.h $(SRCDIR)/graph.h $(SRCDIR)/io.h $(SRCDIR)/align.h $(SRCDIR)/utils.h $(SRCDIR)/trie.h $(EXTDIR)/GraphAligner/GfaGraph.h
DEPS = $(patsubst %, %, $(_DEPS))

_OBJ = $(SRCDIR)/argparse.o $(SRCDIR)/astar.o $(SRCDIR)/gfa2graph.o $(SRCDIR)/align.o $(SRCDIR)/trie.o $(EXTDIR)/GraphAligner/GfaGraph.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

$(shell mkdir -p bin)
$(shell mkdir -p obj/src)
$(shell mkdir -p obj/ext/GraphAligner)

$(ASTARIX_BIN): $(SRCDIR)/astarix.cpp $(OBJ)
	$(GPP) $< -o $@ $(OBJ) $(LINKFLAGS)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

test: $(ASTARIX_BIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIX_BIN) align-optimal -t 1 -g $(DATADIR)/ecoli_head10000_linear.gfa -q $(DATADIR)/illumina.fq -o $(TMPDIR)/astar-default
	#$(ASTARIX_BIN) align-optimal -g $(DATADIR)/ecoli_head10000_linear.gfa -q $(DATADIR)/illumina.fq -o $(TMPDIR)/astar-default
	#$(ASTARIX_BIN) align-optimal -a dijkstra -g $(DATADIR)/ecoli_head10000_linear.gfa -q $(DATADIR)/illumina.fq -o $(TMPDIR)/dijkstra-default
	#$(ASTARIX_BIN) align-optimal -a astar-prefix -g $(DATADIR)/ecoli_head10000_linear.gfa -q $(DATADIR)/illumina.fq -D 5 -f 0 -d 10 -c 3 -M 0 -e 0 -S 1 -G 1 -o $(TMPDIR)/astar-custom
                    
.PHONY: all clean

clean:
	rm -rf $(ODIR)/*
	rm -rf $(BINDIR)/*
	rm -rf $(TMPDIR)/*
