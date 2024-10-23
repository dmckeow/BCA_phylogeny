# Pipeline in development for Biodiversity Cell Atlas project to perform genome-wide protein orthology inference, sequence clustering, alignment, and phylogeny 

- make the search outputs obvious from the command line name  
- mmseqs2 clustering  
- easy-phylo: implement intermediate file checks   

Installation:
```bash
git clone --recurse-submodules https://github.com/zolotarovgl/phylogeny.git
```

Examples: 
```bash
python main.py search -f data/sample.fasta -g data/genefam.tsv Myosin # creates results_annotation/searches/myo.Myosin.domains.fasta

mkdir -p results_annotation/alignments
python main.py align -f results_annotation/searches/myo.Myosin.domains.fasta -o results_annotation/alignments/test.aln -c 10
mkdir -p results_annotation/gene_trees
python main.py phylogeny -f results_annotation/alignments/test.aln -o results_annotation/gene_trees/test -c 15
```
# Dean new:
## Setup
Clone this pipeline and Broccoli:
```bash
git clone https://github.com/dmckeow/BCA_phylogeny.git
cd BCA_phylogeny
git clone https://github.com/rderelle/Broccoli.git
```
## Pre-pipeline steps (search, cluster - initial preparation of fastas which may or may not be inputs for main pipeline):
To install:
```bash
mamba env create -f envs/search_cluster.yaml
```
Used the same way as before (everything else will be moved into a snakemake pipeline):
```bash
conda activate search_cluster
mkdir -p results_annotation results_annotation/searches results_annotation/clusters
python main.py search -f data/sample.fasta -g data/genefam.tsv Myosin # creates results_annotation/searches/myo.Myosin.domains.fasta
python main.py cluster -f results_annotation/searches/myo.Myosin.domains.fasta -o results_annotation/clusters/myo.Myosin.domains.dmnd.mcl -i 1.5
```

## The pipeline (alignment onwards)
Before trying to run anything pipeline-related, there are two configuration files which you must understand and edit:  
* **config.yaml**  
This is where you set YOUR parameters for snakemake. This includes input files, tool parameters, and the location of certain resources such as databases
* **profiles/slurm/config.yaml**  
Specific configuration file for SLURM. You MUST edit this to have your correct resources, partition, username, etc

```bash
conda activate snakemake
# on SLURM on login node (main pipeline process runs on login, but everything else gets submitted to nodes):
snakemake --profile profiles/slurm --executor slurm # use --forceall to overwrite previous outputs

 ```