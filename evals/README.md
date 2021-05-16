## Evaluation of AStarix

This Snakemake will download, run, and plot the results from all experiments in the paper below:


### Prerequisites

* astarix, GraphAligner (commit 23a0ddf27331c581f12ed480dc51bebd716fd8db), and PaSGAL (commit 50ad80cd0fd85a2cffec08cbf43462949989ea40) must be in your `$PATH`
* Snakemake
* Python 3 (numpy, pandas, matplotlib, Jupyter)

### Execution

Run `snakemake -p --cores all` in this directory. When rerunning the set of experiments, please run `snakemake clean` first (this will not clear the contents of the `raw` folder).

### Sample commands
* A\* seeds heuristic `astarix align-optimal -a seeds -f reads.fq -g graph.gfa -o outdir >output`
* A\* prefix heuristic `astarix align-optimal -a prefix -f reads.fq -g graph.gfa -o outdir >output`
* Dijkstra `astarix align-optimal -f reads.fq -g graph.gfa -a dijkstra -o outdir >output`
* GraphAligner `Aligner -f reads.fq -g graph.gfa >output`
* PaSGAL `PaSGAL -q reads.fq -r graph.vg -m vg -o output -t 1`
