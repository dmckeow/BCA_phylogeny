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
## Search:
To install everything for search:
```bash
conda env create -f envs/search.yaml
```
Search is used the same way as before (everything else will be moved into a snakemake pipeline):
```bash
conda activate phylogeny_search
python main.py search -f data/sample.fasta -g data/genefam.tsv Myosin # creates results_annotation/searches/myo.Myosin.domains.fasta
```