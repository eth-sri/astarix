GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++20 -O3 -Iext/plog/include/ -Iext/GraphAligner/ -Iext/concurrentqueue/ -Iext/parallel_hashmap/ -Wno-unused-parameter -Wno-missing-field-initializers

DEBUG ?= 0
ifeq ($(DEBUG), 1)
    CPPFLAGS += -DDEBUG -g
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

_DEPS = $(SRCDIR)/argparse.h $(SRCDIR)/dijkstra.h $(SRCDIR)/astar-prefix.h $(SRCDIR)/astar-seeds-exact.h $(SRCDIR)/astar-seeds-approx.h $(SRCDIR)/gfa2graph.h $(SRCDIR)/graph.h $(SRCDIR)/io.h $(SRCDIR)/align.h $(SRCDIR)/utils.h $(SRCDIR)/trie.h $(EXTDIR)/GraphAligner/GfaGraph.h
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

test: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))

	# small
	#$(ASTARIXBIN) align-optimal -a astar-seeds-exact -t 1 $(RUNFLAGS) -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-seeds-exact --fixed_trie_depth 1 --astar_seeds_max_errors 0
	$(ASTARIXBIN) align-optimal -a astar-seeds       -t 1 $(RUNFLAGS) -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-seeds-approx --fixed_trie_depth 1 --astar_seeds_max_errors 0
	$(ASTARIXBIN) align-optimal -a astar-seeds       -t 1 $(RUNFLAGS) -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-seeds-approx --fixed_trie_depth 1 --astar_seeds_max_errors 1
	#$(ASTARIXBIN) align-optimal -a astar-prefix -t 6 $(RUNFLAGS) -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/astar-prefix
	#$(ASTARIXBIN) align-optimal -a dijkstra -g $(DATADIR)/ecoli_head10000_linear/graph.gfa -q $(DATADIR)/ecoli_head10000_linear/illumina.fq -o $(TMPDIR)/ecoli_head10000_linear/dijkstra-default

	#python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head10000_linear/astar-prefix/alignments.tsv $(TMPDIR)/ecoli_head10000_linear/dijkstra-default/alignments.tsv
	#python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head10000_linear/astar-seeds-exact/alignments.tsv $(TMPDIR)/ecoli_head10000_linear/dijkstra-default/alignments.tsv
	#python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head10000_linear/astar-seeds-approx/alignments.tsv $(TMPDIR)/ecoli_head10000_linear/dijkstra-default/alignments.tsv

bigtest:
	# 10000 reads
	$(ASTARIXBIN) align-optimal -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default
	$(ASTARIXBIN) align-optimal -a dijkstra -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/dijkstra-default
	python3 $(TESTSDIR)/compare_profilings.py $(TMPDIR)/ecoli_head1000000_linear/astar-default/alignments.tsv $(TMPDIR)/ecoli_head1000000_linear/dijkstra-default/alignments.tsv

eval100: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-seeds $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0
	$(ASTARIXBIN) align-optimal -a astar-prefix -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-prefix $(RUNFLAGS) --fixed_trie_depth 1
                    
eval150: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-prefix -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina150.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default $(RUNFLAGS) --fixed_trie_depth 0
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina150.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0

eval250: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	#$(ASTARIXBIN) align-optimal -a astar-prefix -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina250.fq -o $(TMPDIR)/ecoli_head1000000_linear_eval250/astar-prefix $(RUNFLAGS) --fixed_trie_depth 0
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina250.fq -o $(TMPDIR)/ecoli_head1000000_linear_eval250/astar-seeds $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0 --astar_seeds_max_indels 5 --astar_seeds_backwards_algo dfs_for_linear
	#$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina250.fq -o $(TMPDIR)/ecoli_head1000000_linear_eval250/astar-seeds $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0 --astar_seeds_max_indels 5 --astar_seeds_backwards_algo bfs

	#$(ASTARIXBIN) align-optimal -a astar-seeds -t 6 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/illumina250.fq -o $(TMPDIR)/ecoli_head1000000_linear/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0 --astar_seeds_max_indels 10

# pbsim --data-type CLR --depth 2 --model_qc ../pbsim_profiles/model_qc_clr ecoli.fasta
eval_long_ccs: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	#$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long_500_ccs/ecoli_10.fastq -o $(TMPDIR)/ecoli_head1000000_linear_long_500_ccs/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 1 -G 1 -S 1  --astar_seeds_backwards_algo dfs_for_linear
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long_300_ccs/ecoli_100.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_300_ccs/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0 -G 1 -S 1  --astar_seeds_backwards_algo dfs_for_linear
	#$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long_300_ccs/ecoli_100.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_300_ccs/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0 -G 1 -S 1  --astar_seeds_backwards_algo bfs

#	#$(ASTARIXBIN) align-optimal -a astar-seeds-exact -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0 $(RUNFLAGS) 
#	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long_1000_ccs/ecoli_100.fastq -o $(TMPDIR)/ecoli_head1000000_linear_long_1000_ccs/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0 -G 1 -S 1  --astar_seeds_backwards_algo dfs_for_linear --astar_seeds_len 14 --tree_depth 12 --astar_seeds_max_indels 15

eval_long_clr: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
#	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long_500_clr/ecoli_10.fastq -o $(TMPDIR)/ecoli_head1000000_linear_long_500_clr/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 2 -G 1 -S 1  --astar_seeds_backwards_algo dfs_for_linear --astar_seeds_len 12 --tree_depth 12 --astar_seeds_max_indels 10
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long_500_clr/ecoli_10.fastq -o $(TMPDIR)/ecoli_head1000000_linear_long_500_clr/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 2 -G 1 -S 1  --astar_seeds_backwards_algo dfs_for_linear --astar_seeds_len 14 --tree_depth 12 --astar_seeds_max_indels 20
#	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long_clr/ecoli_10.fastq -o $(TMPDIR)/ecoli_head1000000_linear_long_clr/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 2 -G 1 -S 1  --astar_seeds_backwards_algo dfs_for_linear --astar_seeds_len 14 --tree_depth 12 --astar_seeds_max_indels 20
	#$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long_1000_clr/ecoli_100.fastq -o $(TMPDIR)/ecoli_head1000000_linear_long_1000_clr/astar-default $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 2 -G 1 -S 1  --astar_seeds_backwards_algo dfs_for_linear --astar_seeds_len 14 --tree_depth 12 --astar_seeds_max_indels 25

#pbsim --data-type CLR --depth 2 --model_qc ../pbsim_profiles/model_qc_clr --accuracy-mean 0.90 --length-min 300 --length-max 300 ecoli.fasta
eval_long300_90: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long300_90.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads/astar-seeds_errors_2 $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 3 -G 1 -S 1 --astar_seeds_len 15 --tree_depth 13

eval_long300_05: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long300_05.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads/astar-seeds_errors_2 $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 2 -G 1 -S 1 --astar_seeds_len 15 --tree_depth 13

eval_manual: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/ecoli_head1000000_linear/graph.gfa -q $(DATADIR)/ecoli_head1000000_linear/long_manual.fq -o $(TMPDIR)/ecoli_head1000000_linear_long_reads_manual/astar-seeds_errors_2 $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 2 -G 1 -S 1 

run_mhc: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
#	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g evals/graphs/pasgal-MHC1.gfa -q evals/reads/M1_reads100.fa -o evals/results/MHC1-astarix-seeds-dfs     $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0 --astar_seeds_backwards_algo dfs_for_linear
#	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g evals/graphs/pasgal-MHC1.gfa -q evals/reads/M1_reads100.fa -o evals/results/MHC1-astarix-seeds-bfs     $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0 --astar_seeds_backwards_algo bfs
#	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g evals/graphs/pasgal-MHC1.gfa -q evals/reads/M1_reads100.fa -o evals/results/MHC1-astarix-seeds-complex $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0 --astar_seeds_backwards_algo complex
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g evals/graphs/pasgal-MHC1.gfa -q evals/reads/M1_reads100.fa -o evals/results/MHC1-astarix-seeds-topsort $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0 --astar_seeds_backwards_algo topsort

eval_chr22: $(ASTARIXBIN)
	$(shell mkdir -p $(TMPDIR))
	$(ASTARIXBIN) align-optimal -a astar-seeds -t 1 -g $(DATADIR)/chr22/HG_22_linear.gfa -q $(DATADIR)/chr22/chr22.fq -o $(TMPDIR)/chr22/astar-seeds $(RUNFLAGS) --fixed_trie_depth 1 --astar_seeds_max_errors 0 -G 5 -S 1 --astar_seeds_backwards_algo dfs_for_linear

.PHONY: all clean

clean:
	#rm -rf $(ODIR)/*
	rm -rf debug/
	rm -rf release/
	rm -rf $(TMPDIR)/
