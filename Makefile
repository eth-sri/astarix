GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++2a -Iext/plog/include/ -Iext/GraphAligner/ -Iext/concurrentqueue/ -Iext/parallel_hashmap/ -Wno-unused-parameter -Wno-missing-field-initializers

DEBUG ?= 0
ifeq ($(DEBUG), 1)
    CPPFLAGS += -DDEBUG -g -O0 -fno-inline
	BINDIR=debug
	RUNFLAGS += -v 2
else
    CPPFLAGS += -DNDEBUG -O3 -march=native
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
MINIMAPBIN=/home/pesho/libs/minimap2-2.17_x64-linux/minimap2
VGBIN=vg
LIBS= #-lm -lz 

_DEPS = $(SRCDIR)/argparse.h $(SRCDIR)/dijkstra.h $(SRCDIR)/astar-prefix.h $(SRCDIR)/astar-seeds.h $(SRCDIR)/gfa2graph.h $(SRCDIR)/graph.h $(SRCDIR)/io.h $(SRCDIR)/align.h $(SRCDIR)/utils.h $(SRCDIR)/trie.h $(EXTDIR)/GraphAligner/GfaGraph.h
DEPS = $(patsubst %, %, $(_DEPS))

_OBJ = $(SRCDIR)/argparse.o $(SRCDIR)/astar-prefix.o $(SRCDIR)/gfa2graph.o $(SRCDIR)/graph.o $(SRCDIR)/io.o $(SRCDIR)/align.o $(SRCDIR)/utils.o $(SRCDIR)/trie.o $(EXTDIR)/GraphAligner/GfaGraph.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

$(shell mkdir -p $(BINDIR))
$(shell mkdir -p $(ODIR)/src)
$(shell mkdir -p $(ODIR)/ext/GraphAligner)

$(ASTARIXBIN): $(SRCDIR)/astarix.cpp $(DEPS) $(OBJ)
	$(GPP) $< -o $@ $(OBJ) $(LINKFLAGS)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

gfa2fa: $(SRCDIR)/gfa2fa.cpp $(DEPS) $(OBJ)
	$(GPP) $< -o $@ $(OBJ) $(LINKFLAGS)

test: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))

	# small
	$(ASTARIXBIN) align-optimal -a astar-seeds       -t 1 $(RUNFLAGS) -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-seeds --fixed_trie_depth 1 --seeds_len 25
	#$(ASTARIXBIN) align-optimal -a astar-prefix -t 6 $(RUNFLAGS) -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-prefix
	#$(ASTARIXBIN) align-optimal -a dijkstra -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/dijkstra-default

	#python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head10000_linear/astar-prefix/alignments.tsv $(TMPDIR)/ecoli_head10000_linear/dijkstra-default/alignments.tsv

bigtest:
	# 10000 reads
	$(ASTARIXBIN) align-optimal -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default
	$(ASTARIXBIN) align-optimal -a dijkstra -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/dijkstra-default
	python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head1000000_linear/astar-default/alignments.tsv $(TMPDIR)/ecoli_head1000000_linear/dijkstra-default/alignments.tsv

eval100: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-seeds $(RUNFLAGS) --fixed_trie_depth 1
	$(ASTARIXBIN) align-optimal -a astar-prefix -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-prefix $(RUNFLAGS) --fixed_trie_depth 1
                    
eval150: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-prefix -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina150.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-prefix $(RUNFLAGS) --fixed_trie_depth 0
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina150.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-seeds $(RUNFLAGS) --fixed_trie_depth 1
	
	$(shell mkdir -p $(TMPDIR)/ecoli_head1000000_linear/vg/)
	$(VGBIN) map -d $(DATADIR)/ecoli_head1000000_linear/ecoli_1M -f $(DATADIR)/ecoli_head1000000_linear/illumina150.fq -q 1 -z 1 -o 1 -y 1 -L 0 -j -N 1 >$(TMPDIR)/ecoli_head1000000_linear/vg/illumina150.json 

# art 
eval250: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina_250_10000.fq -o $(TMPDIR)/ecoli_head1000000_linear_eval250/astar-seeds $(RUNFLAGS) --fixed_trie_depth 1
	$(ASTARIXBIN) align-optimal -a astar-prefix -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina250.fq -o $(TMPDIR)/ecoli_head1000000_linear_eval250/astar-prefix $(RUNFLAGS) --fixed_trie_depth 1
	$(ASTARIXBIN) align-optimal -a dijkstra -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina250.fq -o $(TMPDIR)/ecoli_head1000000_linear_eval250/dijkstra $(RUNFLAGS) --fixed_trie_depth 1


# pbsim --data-type CLR --depth 2 --model_qc ../pbsim_profiles/model_qc_clr ecoli.fasta
eval_ccs: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(shell mkdir -p $(TMPDIR)/ecoli_head1000000_linear_ccs/vg/)
	$(shell mkdir -p $(TMPDIR)/ecoli_head1000000_linear_ccs/minimap2/)

	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/ecoli.fasta -q $(DATADIR)/ecoli_head1000000_linear/ccs.fastq -o $(TMPDIR)/ecoli_head1000000_linear_ccs/astar-seeds $(RUNFLAGS) --fixed_trie_depth 1 -G 1 -S 1
	$(MINIMAPBIN) -a -x map-pb $(DATADIR)/ecoli_head1000000_linear/ecoli.fasta $(DATADIR)/ecoli_head1000000_linear/ccs.fastq >$(TMPDIR)/ecoli_head1000000_linear_ccs/minimap2/aln.sam

eval_clr: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(shell mkdir -p $(TMPDIR)/ecoli_head1000000_linear_long_clr/minimap2/)
	$(MINIMAPBIN) -ax map-pb $(DATADIR)/ecoli_head1000000_linear/ecoli.fasta $(DATADIR)/ecoli_head1000000_linear/long_clr/sd_0001.fastq >$(TMPDIR)/ecoli_head1000000_linear_long_clr/minimap2/aln.sam

	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/ecoli.fasta -q $(DATADIR)/ecoli_head1000000_linear/long_500_clr/ecoli_10.fastq -o $(TMPDIR)/ecoli_head1000000_linear_long_500_clr/astar-seeds $(RUNFLAGS) --fixed_trie_depth 1 -G 1 -S 1 --seeds_len 14 --tree_depth 12

#pbsim --data-type CLR --depth 2 --model_qc ../pbsim_profiles/model_qc_clr --accuracy-mean 0.90 --length-min 300 --length-max 300 ecoli.fasta
eval_long300_90: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long300_90.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads/astar-seeds_errors_2 $(RUNFLAGS) --fixed_trie_depth 1 -G 1 -S 1 --seeds_len 15 --tree_depth 13

eval_long300_05: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long300_05.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads/astar-seeds_errors_2 $(RUNFLAGS) --fixed_trie_depth 1 -G 1 -S 1 --seeds_len 15 --tree_depth 13

eval_manual: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long_manual.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads_manual/astar-seeds_errors_2 $(RUNFLAGS) --fixed_trie_depth 1 -G 1 -S 1 

run_mhc: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g evals/graphs/pasgal-MHC1.gfa -q evals/reads/M1_reads100.fa -o evals/results/MHC1-astarix-seeds-topsort $(RUNFLAGS) --fixed_trie_depth 1

chr22_linear: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(shell mkdir -p $(TMPDIR)/chr22_linear/minimap2/)
	$(shell mkdir -p $(TMPDIR)/chr22_linear/vg/)

	$(MINIMAPBIN) -ax sr $(DATADIR)/chr22_linear/chr22.fa $(DATADIR)/chr22_linear/chr22_100.fq >$(TMPDIR)/chr22_linear/minimap2/aln.sam
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/chr22_linear/HG_22_linear.gfa -q $(DATADIR)/chr22_linear/chr22_100.fq -o $(TMPDIR)/chr22_linear/astar-seeds $(RUNFLAGS) --fixed_trie_depth 1 -G 5 -S 1
.PHONY: all clean

clean:
	#rm -rf $(ODIR)/*
	rm -rf debug/
	rm -rf release/
	rm -rf $(TMPDIR)/
